
import numpy as np
from scipy import stats
from numba import njit, prange, set_num_threads
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm


@njit(parallel=True)
def _permute_and_center_batch(X_ranked, batch_size):
    n_samples, n_features = X_ranked.shape
    X_centered = np.empty((batch_size, n_samples, n_features), dtype=np.float32)
    for r in prange(batch_size):
        for j in range(n_features):
            perm = np.random.permutation(X_ranked[:, j])
            col_mean = np.mean(perm)
            for i in range(n_samples):
                X_centered[r, i, j] = perm[i] - col_mean
    return X_centered

def _compute_batch_correlations(X_centered, observed_corr):
    # Vectorized correlation for a batch
    numerators = np.einsum('rij,rik->rjk', X_centered, X_centered)
    stds = np.sqrt(np.einsum('rij,rij->rj', X_centered, X_centered))
    stds[stds == 0] = float('inf')
    denom = np.einsum('ri,rj->rij', stds, stds)
    perms = numerators / denom
    return np.sum(np.abs(perms) >= np.abs(observed_corr), axis=0)

def parallel_spearman_perm_test_low_mem(
    X, 
    n_resamples=1000, 
    batch_size=100, 
    early_stop=False, 
    pval_tol=1e-3, 
    check_every=200,
    threads=4
):
    """
    Memory-efficient, multi-threaded Spearman permutation test with early stopping and progress bar.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Input data matrix.
    n_resamples : int, default=1000
        Total number of permutations.
    batch_size : int, default=100
        Number of permutations to process per batch. Controls memory use
    early_stop : bool, default=False
        If True, enables early stopping when p-values stabilize.
    pval_tol : float, default=1e-3
        Convergence tolerance for maximum change in p-values.
    check_every : int, default=200
        Frequency (in permutations) to check convergence.
    threads : int, default=4
        Number of threads to use for batch correlation computation.

    Returns
    -------
    observed_corr : ndarray of shape (n_features, n_features)
        Observed Spearman correlation matrix.
    pvals : ndarray of shape (n_features, n_features)
        Two-sided p-values based on permutation distribution.
    """
    # print('Ranking samples...', end='\r')
    # t0 = time()
    X_ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=X).astype(np.float32)
    # t1 = time()
    # print(f'Ranked samples in {t1-t0} sec')
    
    # print('Correlating observed features...', end='\r')
    # t0 = time()
    observed_corr = np.corrcoef(X_ranked.T)
    # t1 = time()
    # print(f'Correlated observed features in {t1-t0} sec')

    n_features = X.shape[1]
    exceed_count = np.zeros((n_features, n_features), dtype=np.int32)
    prev_pvals = None
    total_perms = 0

    with ThreadPoolExecutor(max_workers=threads) as executor:
        batch_ranges = range(0, n_resamples, batch_size)
        with tqdm(total=n_resamples, desc="Permuting") as pbar:
            for start in batch_ranges:
                curr_batch = min(batch_size, n_resamples - start)
                X_centered = _permute_and_center_batch(X_ranked, curr_batch)

                # Split batch among threads
                chunk_size = curr_batch // threads
                futures = []
                for t in range(threads):
                    s = t * chunk_size
                    e = (t + 1) * chunk_size if t < threads - 1 else curr_batch
                    chunk = X_centered[s:e]
                    futures.append(executor.submit(_compute_batch_correlations, chunk, observed_corr))

                # Sum contributions from all threads
                for f in futures:
                    exceed_count += f.result()

                total_perms += curr_batch
                pbar.update(curr_batch)

                # Early stopping
                if early_stop and total_perms % check_every == 0:
                    curr_pvals = exceed_count / total_perms
                    if prev_pvals is not None:
                        max_change = np.max(np.abs(curr_pvals - prev_pvals))
                        if max_change < pval_tol:
                            print(f"\nEarly stopping at {total_perms} permutations (max pval change: {max_change:.2e})")
                            break
                    prev_pvals = curr_pvals.copy()

    pvals = exceed_count / total_perms
    np.fill_diagonal(pvals, 0)
    return observed_corr, pvals


def run_benchmark(func, its=10, n_samples=100, n_features=500, n_resamples=1000, **kwargs):
    # Generate random input data
    X = np.random.rand(n_samples, n_features)

    # Commpile for JIT
    func(X, n_resamples=n_resamples, **kwargs)
    
    # Time 
    dts = []
    for i in range(its):
        t0 = time()
        func(X, n_resamples=n_resamples, **kwargs)
        t1 = time()
        dts.append(t1-t0)
        
    return np.mean(dts)



def zeroone(vals,lv,hv,func, **kwargs):
    return( 
        (func(vals, **kwargs) - func(lv, **kwargs)) 
        / (func(hv, **kwargs) - func(lv, **kwargs))
    )

def invfun(x, exp):
    return 1/(x**exp)

def signed_spearman_similarity(r_matrix, p_matrix, method='zscore', exp=1, minp=1e-15, distance=False):
    """
    Compute a signed Spearman similarity matrix using correlation and p-values.

    Parameters
    ----------
    r_matrix : ndarray of shape (n, n)
        Matrix of Spearman correlation coefficients.
    p_matrix : ndarray of shape (n, n)
        Matrix of p-values corresponding to the correlations.
    method : {'zscore', '1-p', 'log','expit', 'inv'}, optional
        Method to weight the correlation 
        (all weights are normalized between 0 (p=1-minp) and 1 (p=minp)):
        - 'zscore': weight = normalized_0_to_1( -Φ⁻¹(p/2) )
        - '1-p': weight = 1 - p
        - 'log': weight = normalized_0_to_1( -log(p) )
        - 'expit': weight = 1/(1+np.exp((p-0.1)*60)) 
            (This is a smooth step from 1 to 0 centered at p=0.1)
        - 'inv' : weight = normalized_0_to_1( 1/(x**exp) )
    exp : float, optional
        Exponential to use with 'inv' option
    minp : float, optional
        Minimum allowed p value. Lower p-values are clipped to minp and p-values higher
        than 1-minp are clipped to 1-minp
    distance : bool, optional
        Whether or not to return a distance matrix rather than a similarity matrix

    Returns
    -------
    sim_matrix : ndarray of shape (n, n)
        Signed similarity matrix in the range [-1, 1], where sign comes from r,
        and magnitude increases with statistical confidence and magnitude of r.
    - or - (if distance=True)
    dist_matrix : ndarray of shape (n, n)
        Distance matrix in the range [0, 2], where low values are positively
        correlated with high confidence, medium values are uncorrelated or 
        low confidence, and high values are negatively correlated with 
        high confidence.
    """

    r = np.copy(r_matrix)
    p = np.copy(p_matrix)

    # Avoid invalid values
    maxp = 1-minp
    p = np.clip(p, minp, maxp)

    if method == 'zscore':
        # Convert p-values to two-tailed z-scores
        weight = zeroone(p, maxp, minp, stats.norm.isf)
    elif method == 'expit':
        weight = 1/(1+np.exp((p-0.1)*60))
    elif method == '1-p':
        weight = 1 - p
    elif method == 'log':
        weight = zeroone(p, maxp, minp, np.log)
    elif method == 'inv':
        zeroone(p, maxp, minp, invfun, exp=exp)
    else:
        raise ValueError("Unknown method: choose from 'zscore', 'expit', '1-p', 'log', or 'inv'.")

    # Signed similarity: preserves sign of r and modulates by significance
    out_matrix = r * weight

    # DIstance matrix 
    if distance:
        out_matrix = 1 - out_matrix
        
    return out_matrix




set_num_threads(20)

@njit(parallel=True)
def _permute_and_center_nf(X_ranked, n_resamples):
    n_samples, n_features = X_ranked.shape
    X_centered = np.empty((n_resamples, n_samples, n_features))
    # Parallel on n_features
    for j in prange(n_features):
        for r in range(n_resamples):
            perm = np.random.permutation(X_ranked[:, j])
            col_mean = np.mean(perm)
            for i in range(n_samples):
                X_centered[r, i, j] = perm[i] - col_mean
    return X_centered

@njit(parallel=True)
def _permute_and_center_nr(X_ranked, n_resamples):
    n_samples, n_features = X_ranked.shape
    X_centered = np.empty((n_resamples, n_samples, n_features))
    # Parallel on n_resamples
    for r in prange(n_resamples):
        for j in range(n_features):
            perm = np.random.permutation(X_ranked[:, j])
            col_mean = np.mean(perm)
            for i in range(n_samples):
                X_centered[r, i, j] = perm[i] - col_mean
    return X_centered

# Chunked correlation computation
def _compute_chunked_correlations(X_centered, chunk_size=100):
    n_resamples = X_centered.shape[0]
    # indices
    chunks = [(i, min(i + chunk_size, n_resamples)) for i in range(0, n_resamples, chunk_size)]
    
    def compute_chunk(start, end):
        chunk = X_centered[start:end]
        # Dot product
        numerators = np.einsum('rij,rik->rjk', chunk, chunk)
        # feature Standard deviation
        stds = np.sqrt(np.einsum('rij,rij->rj', chunk, chunk))
        stds[stds == 0] = float('inf')  # features where std is 0, make correlation zero
        # Multiply standard deviations
        denominators = np.einsum('ri,rj->rij', stds, stds)
        # return correlation 
        return numerators / denominators

    results = []
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(compute_chunk, start, end) for start, end in chunks]
        for future in futures:
            results.append(future.result())

    return np.concatenate(results, axis=0)

# Main function
def parallel_spearman_permutation_test(X, n_resamples=1000, permute_parallel_on='n_resamples', corr_chunk_size=100):
    """
    Perform a fast, parallelized permutation test for Spearman correlation using Numba and multi-threaded chunked computation.

    This function computes the Spearman correlation matrix for a dataset and estimates p-values for each pairwise correlation
    via a permutation test. It uses Numba to accelerate permutation and centering of ranked data, and threads to parallelize
    the correlation computation by dividing the permutations into chunks.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Input data matrix where rows are samples and columns are features.

    n_resamples : int, default=1000
        Number of permutations to use when building the null distribution of correlations.

    permute_parallel_on : {'n_resamples', 'n_features'}, default='n_resamples'
        Determines the axis along which to parallelize the permutation process using Numba.
        - 'n_resamples': parallelizes across permutation replicates.
        - 'n_features': parallelizes across features.

    corr_chunk_size : int, default=100
        Size of chunks to divide the permutations into for parallel correlation computation
        using threads. Larger chunk sizes reduce thread overhead but use more memory.

    Returns
    -------
    observed_corr : ndarray of shape (n_features, n_features)
        Observed Spearman correlation matrix computed from the ranked input data.

    pvals : ndarray of shape (n_features, n_features)
        Matrix of two-sided p-values estimated from the permutation test for each feature pair.
        Diagonal entries are set to 0.
    """
    # Which dimension to parallelize permutation on?
    # NOTE: numba may not even be parallelizing this
    permfunc = {
        'n_resamples': _permute_and_center_nr,
        'n_features': _permute_and_center_nf,
    }[permute_parallel_on]

    # For spearman correlation, rank the input, then calculate pearson correlation
    X_ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=X)
    # Numba-accelerated (slightly) permutation and centering
    X_centered = permfunc(X_ranked, n_resamples)

    # Parallelized, chunked correlation computation
    permuted_corrs = _compute_chunked_correlations(X_centered, chunk_size=corr_chunk_size)

    # Get observed rho values
    observed_corr = np.corrcoef(X_ranked.T)

    # Calculate probability of observed values from permuted distribution
    pvals = np.mean(np.abs(permuted_corrs) >= np.abs(observed_corr), axis=0)
    # Assign self perfect correlation
    np.fill_diagonal(pvals, 0)

    return observed_corr, pvals




# Main function using pre-ranked input
def fast_vectorized_numba_spearman_permutation_test(X, n_resamples=1000, parallel_on='n_resamples'):
    """
    Perform a fast, vectorized permutation test for Spearman correlation across features.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Input data matrix.
    n_resamples : int, default=1000
        Number of permutations for the null distribution.
    parallel_on : {'n_resamples', 'n_features'}
        Selects which dimension should be parallelized.

    Returns
    -------
    observed_corr : ndarray of shape (n_features, n_features)
        Observed Spearman correlation matrix.
    pvals : ndarray of shape (n_features, n_features)
        Matrix of permutation-based two-sided p-values.
    """
    # Which dimension to parallelize on?
    permfunc = {
        'n_resamples':_permute_and_center_nr,
        'n_features':_permute_and_center_nf,
    }[parallel_on]

    # # Warm up the JIT
    # dtype = X.dtype
    # dummy = np.random.rand(10, 5).astype(dtype)  # small shape but same dimensionality
    # _ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=dummy)
    # _ = permfunc(_ranked, 2)  # small n_resamples

    # Rank real input 
    X_ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=X)

    # Numba-accelerated permutation and centering
    X_centered = permfunc(X_ranked, n_resamples)

    # Vectorized computation of correlations
    numerators = np.einsum('rij,rik->rjk', X_centered, X_centered)
    stds = np.sqrt(np.einsum('rij,rij->rj', X_centered, X_centered))
    stds[stds == 0] = float('inf')  # features where std is 0, make correlation zero
    denominators = np.einsum('ri,rj->rij', stds, stds)
    permuted_corrs = numerators / denominators

    observed_corr = np.corrcoef(X_ranked.T)
    pvals = np.mean(np.abs(permuted_corrs) >= np.abs(observed_corr), axis=0)
    np.fill_diagonal(pvals, 0)

    return observed_corr, pvals


def fast_vectorized_spearman_permutation_test(X, n_resamples=1000):
    """
    Perform a fast, vectorized permutation test for Spearman correlation across features.

    This function computes the observed Spearman correlation matrix between all pairs of features 
    in the input data matrix `X`, and estimates p-values by comparing the observed correlations 
    against a distribution of correlations obtained by permuting the ranks of each feature independently.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Input data matrix where each column represents a variable (feature) and each row is an observation.
    n_resamples : int, default=1000
        Number of random permutations to perform for the null distribution.

    Returns
    -------
    observed_corr : ndarray of shape (n_features, n_features)
        The observed Spearman rank correlation matrix between features.
    pvals : ndarray of shape (n_features, n_features)
        The matrix of two-sided p-values estimated via permutation testing. Each entry represents
        the probability of observing a correlation as or more extreme than the actual value under the null.
        Diagonal elements are set to 0.
    
    Notes
    -----
    - This implementation vectorizes permutation generation and correlation computation for speed.
    - The null distribution is generated by independently permuting the ranks of each feature.
    - Correlations are computed using dot products of centered rank vectors (equivalent to Pearson on ranks).

    Examples
    --------
    observed_corr, pvals = fast_vectorized_spearman_permutation_test(X, n_resamples=1000)
    """
    n_samples, n_features = X.shape

    # Pre-rank
    X_ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=X)  # shape: (n_samples, n_features)

    # Shape: (n_resamples, n_samples, n_features)
    X_permuted = np.stack([
        np.column_stack([np.random.permutation(X_ranked[:, j]) for j in range(n_features)])
        for _ in range(n_resamples)
    ]) 

    # Center each sample (along sample axis)
    X_centered = X_permuted - X_permuted.mean(axis=1, keepdims=True)

    # Compute dot products for numerator Shape is (n_resamples, n_features, n_features)
    numerators = np.einsum('rij,rik->rjk', X_centered, X_centered)

    # Compute standard deviations (diagonal of covariance matrix)
    stds = np.sqrt(np.einsum('rij,rij->rj', X_centered, X_centered))  # shape: (n_resamples, n_features)

    # Denominator outer product per resample
    denominators = np.einsum('ri,rj->rij', stds, stds)  # shape: (n_resamples, n_features, n_features)

    # Final correlation matrix for all permutations
    permuted_corrs = numerators / denominators  # shape: (n_resamples, n_features, n_features)

    # Observed 
    observed_corr = np.corrcoef(X_ranked.T)

    pvals = np.mean(np.abs(permuted_corrs) >= np.abs(observed_corr), axis=0)
    np.fill_diagonal(pvals, 0)  # or 1, depending on interpretation

    return observed_corr, pvals

def fast_spearman_permutation_test(X, n_resamples=1000, random_state=None):
    rng = np.random.default_rng(random_state)
    n_samples, n_features = X.shape

    # Pre-rank the data
    X_ranked = np.apply_along_axis(stats.rankdata, axis=0, arr=X)
    observed = np.corrcoef(X_ranked.T)
    perms = np.stack([
        np.column_stack([np.random.permutation(X_ranked[:, j]) for j in range(n_features)])
        for _ in range(n_resamples)
    ])
    # Allocate memory for permutations
    null_distributions = np.zeros((n_resamples, n_features, n_features))

    for i in range(n_resamples):
        X_ranked_perm = np.column_stack(
            [np.random.permutation(X_ranked[:, j]) for j in range(n_features)]
        )
        corr = np.corrcoef(X_ranked_perm.T)
        null_distributions[i] = corr

    # Compute two-sided p-values
    pvals = np.mean(np.abs(null_distributions) >= np.abs(observed), axis=0)
    return observed, pvals


