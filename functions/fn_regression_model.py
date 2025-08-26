import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder, PolynomialFeatures
import statsmodels.api as sm
from joblib import Parallel, delayed
import warnings
import time
import statsmodels.stats.multitest as smm
from collections import defaultdict


def adjust_pvalues_for_multiple_testing(results, method='fdr_bh'):
    # Gather all p-values and their keys
    all_pvals = []
    keys = []
    for response_col, result in results.items():
        pval_dict = result['perm_pvalues'].get('p_value', {})
        for var, pval in pval_dict.items():
            all_pvals.append(pval)
            keys.append((response_col, var))
    if not all_pvals:
        return results  # No p-values to adjust

    # Adjust p-values
    adjusted_pvals = smm.multipletests(np.array(all_pvals), method=method)[1]

    # Insert adjusted p-values back into results
    for (response_col, var), adj_pval in zip(keys, adjusted_pvals):
        if 'adj_perm_pvalues' not in results[response_col]:
            results[response_col]['adj_perm_pvalues'] = {'p_value': {}}
        results[response_col]['adj_perm_pvalues']['p_value'][var] = adj_pval
    return results


def likelihood_ratio_statistic(y, X_poly_perm, X_onehot, family):
    """
    Compute the likelihood ratio statistic (difference in deviance)
    between a full GLM and a null (non-perm) GLM.
    
    Parameters:
        y: array-like, response variable
        X_full: array-like, full design matrix (including intercept)
        family: statsmodels family object (e.g., sm.families.Gaussian())
    
    Returns:
        lr_stat: float, likelihood ratio statistic (null deviance - full deviance)
    """
    # recompose design matrix
    X_perm = pd.concat([X_onehot, X_poly_perm], axis=1)
    # Fit the full model
    model_full = sm.GLM(y, X_perm, family=family).fit()
    deviance_full = model_full.deviance

    # Fit the null model (intercept only)
    # X_null = sm.add_constant(np.ones(len(y)))
    model_null = sm.GLM(y, X_onehot, family=family).fit()
    deviance_null = model_null.deviance

    # Likelihood ratio statistic
    lr_stat = deviance_null - deviance_full
    return model_full, lr_stat



def permute_between_groups_of_replicate_samples(series, group_df, rng):
    # Combine series and group_df
    df = pd.concat([series.rename('perm_val'), group_df], axis=1)
    group_cols = group_df.columns.tolist()
    # Assign a unique group id
    df['_group_id'] = df.groupby(group_cols, sort=False).ngroup()
    # Get one value per group
    unique_vals = df.groupby('_group_id', sort=False)['perm_val'].first().values
    # Shuffle
    rng.shuffle(unique_vals)
    # Map shuffled values back to groups
    group_to_val = dict(zip(df['_group_id'].unique(), unique_vals))
    df['perm_val'] = df['_group_id'].map(group_to_val)
    return df['perm_val'].rename(series.name)


def permute_within_groups_02(series, group_df, rng):
    df = pd.concat([series, group_df], axis=1)
    group_cols = group_df.columns.tolist()
    permuted = []
    for _, group in df.groupby(group_cols, sort=False):
        permuted_vals = permute_between_groups_of_replicate_samples(
            group[series.name], group, rng
        )
        permuted.append(pd.Series(permuted_vals))
    return pd.concat(permuted)

def remove_groupwise_outliers(y, X, onehot_feature_names, iqr_coefficient=1.5):
    y = pd.Series(y, index=X.index)
    if onehot_feature_names:
        df = pd.concat([y.rename('y'), X], axis=1)
        keep_idx = []
        for _, group in df.groupby(onehot_feature_names, sort=False):
            q1 = group['y'].quantile(0.25)
            q3 = group['y'].quantile(0.75)
            iqr = q3 - q1
            lower = q1 - iqr_coefficient * iqr
            upper = q3 + iqr_coefficient * iqr
            mask = (group['y'] >= lower) & (group['y'] <= upper)
            keep_idx.extend(group.index[mask])
        keep_idx = sorted(set(keep_idx))
    else:
        q1 = y.quantile(0.25)
        q3 = y.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - iqr_coefficient * iqr
        upper = q3 + iqr_coefficient * iqr
        mask = (y >= lower) & (y <= upper)
        keep_idx = sorted(set(X.index[mask]))
    n_outliers = y.shape[0] - len(keep_idx)
    return keep_idx, n_outliers

def preprocess_onehot_matrix(data, covariate_cols, poly_cols=None):
    """Precompute one-hot encoding for categorical covariates only."""
    if poly_cols is None:
        poly_cols = []
    categorical_cols = [col for col in covariate_cols if col not in poly_cols]
    if not categorical_cols:
        return pd.DataFrame(index=data.index), [], []
    ohe = OneHotEncoder(drop='first', sparse_output=False)
    X_onehot = ohe.fit_transform(data[categorical_cols])
    feature_names = ohe.get_feature_names_out(categorical_cols)
    return pd.DataFrame(X_onehot, columns=feature_names, index=data.index), list(feature_names), categorical_cols

def compute_poly_matrix(data, poly_cols, poly_degree=2):
    """Compute polynomial features for continuous covariates."""
    if not poly_cols:
        return pd.DataFrame(index=data.index), [], {}
    poly = PolynomialFeatures(degree=poly_degree, include_bias=False)
    X_poly = poly.fit_transform(data[poly_cols])
    poly_feature_names = poly.get_feature_names_out(poly_cols)
    # Build poly_groups: for each base var, collect its derived features
    poly_groups = {}
    for base_var in poly_cols:
        base_var_features = [fname for fname in poly_feature_names if fname.startswith(base_var)]
        poly_groups[base_var] = base_var_features
    return pd.DataFrame(X_poly, columns=poly_feature_names, index=data.index), list(poly_feature_names), poly_groups

def estimate_alpha(y, X, thresh_disp=1.5):
    X = sm.add_constant(X)
    poisson_model = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    mu = poisson_model.fittedvalues
    mu = np.where(mu <= 0, 1e-8, mu)
    # Calculate dispersion ratio
    pearson_chi2 = np.sum((y - mu) ** 2 / mu)
    df_resid = poisson_model.df_resid
    dispersion_ratio = pearson_chi2 / df_resid
    if dispersion_ratio > thresh_disp:
        resid_pearson_sq = ((y - mu) ** 2 - y) / mu
        mask = np.isfinite(resid_pearson_sq) & np.isfinite(mu)
        if not np.any(mask):
            raise ValueError("No finite values to estimate alpha.")
        aux_olsr = sm.OLS(resid_pearson_sq[mask], mu[mask]).fit()
        alpha = aux_olsr.params[0]
    else:
        alpha = 0
    return alpha, pearson_chi2, dispersion_ratio


# def fit_single_model_glm_nb_permute_predictor_grouped(
#     X_onehot, X_poly, data, y, alpha, perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, n_perm, stat, random_state
# ):
#     timings = {}
#     timings['permutation_steps'] = {}  # Will store timings for each variable

#     X_onehot = sm.add_constant(X_onehot)
#     group_df = X_onehot[onehot_feature_names] if onehot_feature_names else None
#     X = pd.concat([X_onehot, X_poly], axis=1)
#     model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
#     perm_pvalues = pd.DataFrame()
#     perm_dict = {}
#     n_skipped = 0

#     if perm_test_vars:
#         for var in perm_test_vars:
#             # Determine if var is a polynomial feature
#             if poly_cols:
#                 for base in poly_cols:
#                     if var in poly_groups.get(base, []):
#                         var = base
#                         break
#             coef_obs = model.params[var] if stat == 'coef' else model.tvalues[var]
#             perm_stats = []
#             # Timing for this variable's permutations
#             step_timings = {
#                 'permute_data': [],
#                 'compute_poly': [],
#                 'recompose_X': [],
#                 'fit_model': [],
#                 'total': []
#             }
#             for i in range(n_perm):
#                 t_total = time.time()
#                 # 1. Permute data
#                 t0 = time.time()
#                 rng = np.random.default_rng(random_state)
#                 data_perm = data.copy()
#                 if group_df is not None and not group_df.empty:
#                     data_perm[var] = permute_within_groups_02(data_perm[var], group_df, rng)
#                 else:
#                     data_perm[var] = permute_between_groups_of_replicate_samples(
#                         data_perm[var], data_perm, rng
#                     )
#                 step_timings['permute_data'].append(time.time() - t0)

#                 # 2. Recompute polynomial features
#                 t0 = time.time()
#                 X_poly_perm, _, _ = compute_poly_matrix(data_perm, poly_cols, poly_degree)
#                 step_timings['compute_poly'].append(time.time() - t0)

#                 # 3. Recompose design matrix
#                 t0 = time.time()
#                 X_perm = pd.concat([X_onehot, X_poly_perm], axis=1)
#                 X_perm = sm.add_constant(X_perm, has_constant='add')
#                 step_timings['recompose_X'].append(time.time() - t0)

#                 # 4. Fit permuted model
#                 t0 = time.time()
#                 try:
#                     model_perm = sm.GLM(y, X_perm, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
#                     perm_stat = model_perm.params[var] if stat == 'coef' else model_perm.tvalues[var]
#                     perm_stats.append(perm_stat)
#                 except Exception:
#                     n_skipped += 1
#                     continue
#                 step_timings['fit_model'].append(time.time() - t0)
#                 # 5. Total time for this permutation
#                 step_timings['total'].append(time.time() - t_total)

#             # Store timings for this variable
#             timings['permutation_steps'][var] = step_timings

#             p_value = (np.abs(perm_stats) >= np.abs(coef_obs)).mean() if perm_stats else np.nan
#             perm_pvalues.loc[var, 'p_value'] = p_value
#             perm_dict[var] = perm_stats

#     return model, perm_pvalues, perm_dict, alpha, n_skipped, timings


def fit_single_model_glm_nb_permute_predictor_grouped(
    X_onehot, X_poly, data, y, alpha, perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, n_perm, stat, random_state
):


    X = pd.concat([X_onehot, X_poly], axis=1)
    # Get model family
    if alpha > 0:
        family = sm.families.NegativeBinomial(alpha=alpha) 
    else:
        family = sm.families.Poisson()
    # Fit model
    model = sm.GLM(y, X, family=family).fit()
    # model, stat = likelihood_ratio_statistic(y, X_poly, X_onehot, family)
    # fit the null model (no poly)
    model_null = sm.GLM(y, X_onehot, family=family).fit()
    deviance_null = model_null.deviance
    lr_stat = deviance_null - model.deviance

    # Set up permutation outputs
    perm_pvalues = pd.DataFrame()
    perm_dict = defaultdict(dict)
    rng = np.random.default_rng(random_state)
    n_skipped = 0
    # Set up groups to permute within
    group_df = X_onehot[onehot_feature_names] if onehot_feature_names else None
    # Permute variables to get coefficient distribution
    if perm_test_vars:
        for var in perm_test_vars:
            # Determine if var is a polynomial feature
            if poly_cols:
                for base in poly_cols:
                    if var in poly_groups.get(base, []):
                        var = base
                        break
            # coef_obs = model.params[var] if stat == 'coef' else model.tvalues[var]
            perm_stats = defaultdict(list)
            for i in range(n_perm):
                # Permute base variable (for polynomials) or var (for non-poly)
                data_perm = data.copy()
                if group_df is not None and not group_df.empty:
                    data_perm[var] = permute_within_groups_02(data_perm[var], group_df, rng)
                else:
                    # data_perm[var] = rng.permutation(data_perm[var].values)
                    data_perm[var] = permute_between_groups_of_replicate_samples(
                        data_perm[var], data_perm, rng
                    )
                X_poly_perm, _, _ = compute_poly_matrix(data_perm, poly_cols, poly_degree)
                # # Recompose design matrix
                X_perm = pd.concat([X_onehot, X_poly_perm], axis=1)
                # X_perm = sm.add_constant(X_perm, has_constant='add')
                try:
                    model_perm = sm.GLM(y, X_perm, family=family).fit()
                    perm_stat = deviance_null - model_perm.deviance
                        # perm_stat = model_perm.params[var] if stat == 'coef' else model_perm.tvalues[var]
                    # _, perm_stat = likelihood_ratio_statistic(y, X_poly_perm, X_onehot, family)
                    perm_stats['lr_stat'].append(perm_stat)
                    for c in X_poly_perm.columns:
                        perm_stats[c].append(model_perm.params[c])
                except Exception:
                    n_skipped += 1
                    continue

            # define main model stat values
            dict_statname_modelval = {'lr_stat': lr_stat}
            for c in X_poly.columns:
                dict_statname_modelval[c] = model.params[c]
            # Compute permutation p values
            for statname, permvals in perm_stats.items():
                modelval = dict_statname_modelval[statname]
                p_value = (np.abs(permvals) >= np.abs(modelval)).mean() if permvals else np.nan
                perm_pvalues.loc[statname, 'p_value'] = p_value
                perm_dict[var][statname] = permvals

    return model, lr_stat, perm_pvalues, perm_dict, alpha, n_skipped


# def process_response_col(
#     response_col, X_onehot, X_poly, data, perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, n_perm, stat, random_state
# ):
#     timings = {}
#     with warnings.catch_warnings(record=True) as wlist:
#         warnings.simplefilter('always')

#         # Time: Response values extraction
#         t0 = time.time()
#         y = data[response_col].values.copy()
#         timings['extract_response'] = time.time() - t0

#         # Time: Remove outliers
#         t0 = time.time()
#         keep_idx = remove_groupwise_outliers(y, X_onehot, onehot_feature_names)
#         y = pd.Series(y, index=X_onehot.index).loc[keep_idx].values
#         X_onehot = X_onehot.copy().loc[keep_idx]
#         X_poly = X_poly.copy().loc[keep_idx]
#         data_sub = data.copy().loc[keep_idx]
#         timings['remove_outliers'] = time.time() - t0

#         # Time: Estimate alpha
#         t0 = time.time()
#         X = pd.concat([X_onehot, X_poly], axis=1)
#         alpha = estimate_alpha(y, X)
#         timings['estimate_alpha'] = time.time() - t0

#         # Time: Fit model and permutations
#         t0 = time.time()
#         model, perm_pvalues, perm_dict, alpha_used, n_skipped, timings_glm = fit_single_model_glm_nb_permute_predictor_grouped(
#             X_onehot.copy(), X_poly.copy(), data_sub, y, alpha, 
#             perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, 
#             n_perm, stat, random_state
#         )
#         timings['fit_model_and_permute'] = time.time() - t0

#         return {
#             'result': [response_col, model, perm_pvalues, perm_dict, alpha_used, n_skipped],
#             'error': None,
#             'warnings': [str(w.message) for w in wlist],
#             'timings': {**timings, **timings_glm}
#         }


def process_response_col(
    response_col, X_onehot, X_poly, data, 
    perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, 
    n_perm, stat, random_state, iqr_coefficient
):
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter('always')

        # Response values
        y = data[response_col].values.copy()

        # Remove outliers for each group defined by the onehot encoding
        keep_idx, n_outliers = remove_groupwise_outliers(
            y, X_onehot, onehot_feature_names, iqr_coefficient=iqr_coefficient
        )
        y = pd.Series(y, index=X_onehot.index).loc[keep_idx].values
        X_onehot = X_onehot.copy().loc[keep_idx]
        X_poly = X_poly.copy().loc[keep_idx]
        data_sub = data.copy().loc[keep_idx]

        try:
        # Get alpha for this taxon
            X = pd.concat([X_onehot, X_poly], axis=1)

            alpha, pearson_chi2, dispersion_ratio = estimate_alpha(y, X)
            # Fit taxon
            model, lr_stat, perm_pvalues, perm_dict, alpha_used, n_skipped = fit_single_model_glm_nb_permute_predictor_grouped(
                X_onehot.copy(), X_poly.copy(), data_sub, y, alpha, 
                perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, 
                n_perm, stat, random_state
            )
            model_out = model.params.to_dict()
            perm_pvalues = perm_pvalues.to_dict()
            return {
                'result': [response_col, model_out, lr_stat, 
                            perm_pvalues, perm_dict, 
                            alpha_used, pearson_chi2, dispersion_ratio, 
                            n_skipped, n_outliers],
                'error': None,
                'warnings': [str(w.message) for w in wlist]
            }
        except Exception as e:
            return {
                'result': None, 
                'error': str(e), 
                'response_col': response_col, 
                'warnings': [str(w.message) for w in wlist]
            }

def fit_many_models(
    data, response_cols, covariate_cols, poly_cols=None, poly_degree=2, 
    iqr_coefficient=1.5, p_value_correction_method='fdr_bh', 
    perm_test_vars=None, n_perm=1000, stat='coef', random_state=None, n_jobs=1
):
    """
    Precompute one-hot encoding, then for each permutation, recompute polynomial features after permuting base variables.
    """
    print('chack')
    # Compute onehot encoding
    X_onehot, onehot_feature_names, _ = preprocess_onehot_matrix(data, covariate_cols, poly_cols)
    X_onehot = sm.add_constant(X_onehot)
    # Compute initial polynomial features
    X_poly, poly_feature_names, poly_groups = compute_poly_matrix(data, poly_cols, poly_degree)
    results = {}
    if n_jobs == 1:
        for response_col in response_cols:
            res = process_response_col(
                response_col, X_onehot, X_poly, data, 
                perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups, 
                n_perm, stat, random_state, iqr_coefficient
            )
            if res['error']:
                print(f"Error for response column {res['response_col']}: {res['error']}")
            else:
                response_col, model, lr_stat, perm_pvalues, perm_dict, alpha_used, pearson_chi2, dispersion_ratio, n_skipped, n_outliers = res['result']
                alpha_info = {
                    'alpha_used':alpha_used, 
                    'pearson_chi2':pearson_chi2, 
                    'dispersion_ratio':dispersion_ratio
                }
                results[response_col] = {
                    'model': model,
                    'lr_stat': lr_stat,
                    'perm_pvalues': perm_pvalues,
                    'perm_dict': perm_dict,
                    'alpha': alpha_info,
                    'n_outliers': n_outliers,
                    'n_skipped': n_skipped,
                    'warnings': res['warnings'],
                }
                    # 'timings': res['timings']
    else:
        tasks = (
            delayed(process_response_col)(
                response_col, X_onehot, X_poly, data, 
                perm_test_vars, onehot_feature_names, poly_cols, poly_degree, poly_groups,
                n_perm, stat, random_state, iqr_coefficient
            )
            for response_col in response_cols
        )
        parallel_results = Parallel(n_jobs=n_jobs, verbose=10)(tasks)
        for i, res in enumerate(parallel_results):
            if res['error']:
                print(f"Error for response column {res['response_col']}: {res['error']}")
            else:
                response_col, model, lr_stat, perm_pvalues, perm_dict, alpha_used, pearson_chi2, dispersion_ratio, n_skipped, n_outliers = res['result']
                alpha_info = {
                    'alpha_used':alpha_used, 
                    'pearson_chi2':pearson_chi2, 
                    'dispersion_ratio':dispersion_ratio
                }
                results[response_col] = {
                    'model': model,
                    'lr_stat': lr_stat,
                    'perm_pvalues': perm_pvalues,
                    'perm_dict': perm_dict,
                    'alpha': alpha_info,
                    'n_outliers': n_outliers,
                    'n_skipped': n_skipped,
                    'warnings': res['warnings']
                }
    # Multiple testing correction
    if p_value_correction_method:
        results = adjust_pvalues_for_multiple_testing(
            results, method=p_value_correction_method
        )
    return results

# Utility functions (permute_within_groups and remove_groupwise_outliers) remain unchanged from your code.
