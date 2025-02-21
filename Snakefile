#!/usr/bin/env -S snakemake -s

"""
Pipeline to extract data from the Armbrust Gradients metaT dataset

"""

# =============================================================================
# Imports
# =============================================================================

import os
import csv
import json
import math
import glob
import tarfile
import fnmatch
import subprocess
import pandas as pd
from collections import defaultdict
os.environ['XDG_DATA_HOME'] = config['XDG_DATA_HOME'] or os.environ['PWD']  # Dir to write the ncbitaxa database (written into 'ete' added to path)
from ete4 import NCBITaxa

## TEMP
sys.path.append('/scratch/bgrodner/repo-armbrust-metat-search')
from functions.cl_tree_trim_02 import TreeTrim


# =============================================================================
# Functions
# =============================================================================

def get_tar_names(fn_tar, re_tar, search_output_dir, fn_metat):
    # Make tar name dir
    out_dir = config['snakemake_output_dir'] + f'/{search_output_dir}/tar_names'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Make filename
    fn_tar_names = os.path.split(fn_tar)[1]
    fn_tar_names = re.sub('.tar.gz','-tar_names.txt', fn_tar_names)
    fn_tar_names_full = f'{out_dir}/{fn_metat}-{fn_tar_names}'
    # Check if the file already exists
    if not os.path.exists(fn_tar_names_full):
        print('Getting tar sub-filenames')
        tar_names = []
        # Read the tar names
        with tarfile.open(fn_tar, "r:gz") as tar:
            for name in tar:
                n = name.name
                if bool(re.search(re_tar, n)):
                    tar_names.append(n)
        # Write file
        with open(fn_tar_names_full, 'w') as f:
            for line in tar_names:
                f.write(f"{line}\n")
    # If it already exists, read the names
    else:
        with open(fn_tar_names_full, 'r') as f:
            tar_names = f.read().splitlines()
    return tar_names

def get_glob_fns(dir_glob, fn_glob):
    if dir_glob:
        return(glob.glob(f'{dir_glob}/{fn_glob}'))
    else:
        return(glob.glob(fn_glob))



def get_fns(fmt, col_dir='search_output_dir', col_fn_metat='fn_metat', col_fn_targets='fn_targets', col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        fn = fmt.format(
            search_output_dir=row[col_dir],
            fn_metat=row[col_fn_metat], 
            fn_targets=row[col_fn_targets],
        )
        if not col_check:
            # Build and append filename
            fns.append(fn)
        else:
            if row[col_check]:
                fns.append(fn)

    return fns

def get_fns_sample(fmt, col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        # Column to check if appending row filenames
        appnd = False
        if not col_check:
            # Build and append filename
            appnd = True
        else:
            if row[col_check]:
                appnd = True
        if appnd:
            # Use the glob expression to get a list of kallisto filenames
            fns_kallisto = get_glob_fns(
                row.dir_kallisto, row.glob_kallisto
            )
            for fnk_ in fns_kallisto:
                # Check if the file is a tarball
                if row.re_kallisto_tar:
                    # Extract the filenames in the tarball 
                    tar_names = get_tar_names(
                        fnk_, 
                        row.re_kallisto_tar, 
                        row.search_output_dir,
                        row.fn_metat
                    )
                    for fnkt_ in tar_names:  
                        fnkt = os.path.split(fnkt_)[1]
                        fnk = os.path.split(fnk_)[1]
                        fnk_merge = f'{fnk}.{fnkt}'
                        # Build and append output filename
                        fns.append(fmt.format(
                            search_output_dir=row.search_output_dir,
                            fn_metat=row.fn_metat, 
                            fn_targets=row.fn_targets,
                            fn_kallisto=fnk_merge
                        ))
                else: 
                    fnk = os.path.split(fnk_)[1]
                    # Build and append filename
                    fns.append(fmt.format(
                        search_output_dir=row.search_output_dir,
                        fn_metat=row.fn_metat, 
                        fn_targets=row.fn_targets,
                        fn_kallisto=fnk
                    ))
    return fns

def get_fns_row_sample(
    fmt, 
    fn_metat, fn_targets, input_table, 
    col_dir='search_output_dir', 
    col_fn_metat='fn_metat', 
    col_fn_targets='fn_targets', 
    col_check=''
):
    fns = []
    row = get_input_table_row(fn_metat, fn_targets, input_table)
    # Column to check if appending row filenames
    appnd = False
    if not col_check:
        # Build and append filename
        appnd = True
    else:
        if row[col_check]:
            appnd = True
    if appnd:
        # Use the glob expression to get a list of kallisto filenames
        fns_kallisto = get_glob_fns(
            row.dir_kallisto, row.glob_kallisto
        )
        for fnk_ in fns_kallisto:
            # Check if the file is a tarball
            if row.re_kallisto_tar:
                # Extract the filenames in the tarball 
                tar_names = get_tar_names(
                    fnk_, 
                    row.re_kallisto_tar, 
                    row.search_output_dir,
                    row.fn_metat
                )
                for fnkt_ in tar_names:  
                    fnkt = os.path.split(fnkt_)[1]
                    fnk = os.path.split(fnk_)[1]
                    fnk_merge = f'{fnk}.{fnkt}'
                    # Build and append output filename
                    fns.append(fmt.format(
                        search_output_dir=row.search_output_dir,
                        fn_metat=row.fn_metat, 
                        fn_targets=row.fn_targets,
                        fn_kallisto=fnk_merge
                    ))
            else: 
                fnk = os.path.split(fnk_)[1]
                # Build and append filename
                fns.append(fmt.format(
                    search_output_dir=row.search_output_dir,
                    fn_metat=row.fn_metat, 
                    fn_targets=row.fn_targets,
                    fn_kallisto=fnk
                ))
    return fns

def get_fns_target(fmt, col_dir='search_output_dir', col_fn_targets='fn_targets', col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        fn = fmt.format(
            search_output_dir=row[col_dir],
            fn_targets=row[col_fn_targets],
        )
        if not col_check:
            # Build and append filename
            fns.append(fn)
        else:
            if row[col_check]:
                fns.append(fn)
    return set(fns)


def get_input_table_value(fn_metat, fn_targets, input_table, column):
    return input_table.loc[
        (
            (input_table.fn_metat == fn_metat) 
            & (input_table.fn_targets == fn_targets)
        ), 
        column
    ].values[0]

def get_input_table_row(fn_metat, fn_targets, input_table):
    return input_table.loc[
        (
            (input_table.fn_metat == fn_metat) 
            & (input_table.fn_targets == fn_targets)
        )
    ].squeeze()



def get_fn_full(fn_metat, fn_targets, input_table, col_dir, col_fn):
    path = get_input_table_value(fn_metat, fn_targets, input_table, col_dir)
    fn = get_input_table_value(fn_metat, fn_targets, input_table, col_fn)
    if path:
        return path + f'/{fn}'
    else:
        return fn

def get_fn_kallisto_full(fn_metat, fn_targets, fn_kallisto, input_table):
    row = get_input_table_row(
        fn_metat, fn_targets, input_table,
    )
    # Split wildcard to get kallisto tarball filename
    if row.re_kallisto_tar:
        regex = fnmatch.translate(row.glob_kallisto)
        regex = regex.rstrip('\Z')
        fn_kallisto = re.match(regex, fn_kallisto)[0]
    # Assemble full filename
    fns_glob = get_glob_fns(
        row.dir_kallisto, row.glob_kallisto
    )
    for fnk_full in fns_glob:
        if fn_kallisto in fnk_full:
            return fnk_full

def get_fn_kallisto_tar(fn_metat, fn_targets, fn_kallisto, input_table):
    row = get_input_table_row(
        fn_metat, fn_targets, input_table,
    )
    # Split wildcard to get kallisto tarball filename
    if row.re_kallisto_tar:
        regex = fnmatch.translate(row.glob_kallisto + '.')
        regex = regex.rstrip('\Z')
        fn_kallisto_tar = re.split(regex, fn_kallisto)[1]
        return fn_kallisto_tar
    # Assemble full filename
    else:
        return ''

# Parse the kallisto sample name (the sub tar name in tarball)
def parse_fn_kallisto_sn(fn, sn_type='', get_columns=False):
    if not get_columns:
        try:
            ass, sample, lat, ammend, timep, depth, size, rep = [''] * 8
            if sn_type == 'G1NS':
                splt = fn.split('.')
                ass, sm_sz, rep = splt[:3]
                sample, sz = sm_sz.split('_',1)
                size = re.sub('_','.',sz)
            elif sn_type == 'G2NS':
                ass, sample, depth, sz, rep, _ = fn.split('.')
                size = re.sub('_','.',sz)
            elif sn_type == 'G3NS':
                meta_fn = input_table.loc[
                    input_table['name'] == 'G3NS', 
                    'fn_sample_metadata'
                ].values[0]
                meta = pd.read_csv(meta_fn)
                sid = os.path.splitext(fn)[0]
                meta_sample = meta.loc[meta['SampleID'] == sid, :]
                depth = str(meta_sample['Depth'].values[0])
                lat = meta_sample['Latitude'].values[0]
                lat = str(round(float(lat),2))
                ass_sm, dp_sz_rep = fn.split('_', 1)
                ass = re.match(r'.+NS', ass_sm)[0]
                sample = re.search(r'UW\d+', ass_sm)[0]
                dp1, dp2, sz, rep, _ = dp_sz_rep.split('.')
                depth = f'{dp1}.{dp2}'
                size = re.sub('_','.',sz)
            elif sn_type == 'G5':
                ass, sample, ammend, timep, rep, _ = fn.split('.')
            elif sn_type == 'D1':
                ass, sm_rep_tp, _, _ = fn.split('.')
                sample, rep, timep = sm_rep_tp.split('_')
            elif sn_type == 'G1PA':
                ass, fn_ = fn.split('.', 1)
                sample, fn_ = fn_.split('_', 1)
                size = re.search(r'.+um', fn_)[0]
                rep = re.search(r'(?<=um)\w+(?=\.)',fn_)[0]
            elif sn_type == 'G2PA':
                _, ass, sample, depth, sz, rep, _, _ = fn.split('.')
                size = re.sub('_','.',sz)
            elif sn_type == 'G3PA.UW':
                meta_fn = input_table.loc[
                    input_table['name'] == 'G3PA.UW', 
                    'fn_sample_metadata'
                ].values[0]
                meta = pd.read_csv(meta_fn)
                ass, sample_, _, _, _, _ = fn.split('.')
                meta_sample = meta.loc[meta['Alias2'] == sample_, :]
                size = str(meta_sample['Filter'].values[0])
                depth = str(meta_sample['Depth'].values[0])
                rep = str(meta_sample['Replicate'].values[0])
                sample_list = str(meta_sample['Alias1'].values[0]).split(' ')
                sample = sample_list[0]
                if '#' in sample_list[1]:
                    sample += sample_list[1]
            elif sn_type == 'G3PA.diel':
                ass1, ass2, sample, rep, _, _, _, _ = fn.split('.')
                ass = f'{ass1}.{ass2}'
            elif sn_type == 'G3PA.PM':
                ass_sm, dp_tp_sz_rp = fn.split('_', 1)
                ass = re.match(r'.+(?=.UW)', ass_sm)[0]
                sample = re.search(r'UW\d+$', ass_sm)[0]
                depth, tp_sz_rp = dp_tp_sz_rp.split('_',1)
                timep, sz_rp = tp_sz_rp.split('.',1)
                size = re.match(r'.+um(?=\.)', sz_rp)[0]
                rep = re.search(r'(?<=um\.)\w+', sz_rp)[0]
            else:
                raise ValueError(
                    f"""
                    Sample name parse type not provided (sn_type_parse_kallisto column in file table)
                    """
                )        
            return [ass, sample, lat, ammend, timep, depth, size, rep]
        except:
            raise ValueError(
                f"""
                Failed to parse filename:
                {fn}
                Using type:
                {sn_type}
                """
            )
    else:
        return ['assembly', 'sample', 'latitude','ammendment', 'timepoint', 'depth', 'size', 'rep']

def get_fn_list_contigs_diamond(fn_metat, fn_targets, input_table, fn_6tr, fn_std):
    type_contig_name = get_input_table_value(
            fn_metat, fn_targets, input_table, 'type_diamond_contig_name'
        )
    if type_contig_name == '6tr':
        return fn_6tr
    else: 
        return fn_std



# def get_fn_metat_full(fn_metat, fn_targets, input_table):
#     path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_metat')
#     if path:
#         return path + f'/{fn_metat}'
#     else:
#         return fn_metat


# def get_fn_targets_full(fn_metat, fn_targets, input_table):
#     path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_targets')
#     if path:
#         return path + f'/{fn_targets}'
#     else:
#         return fn_targets


# =============================================================================
# Setup
# =============================================================================


input_table = pd.read_csv(
    config['input_table'],
    keep_default_na=False
)

# format start files
fmt_start = (
    config['snakemake_output_dir'] + '/start_files/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}.txt'
)
# format get_KO_contig_dict outputs
fmt_dict = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dict.json'
)
# format dict_to_list outputs
fmt_contig_list_6tr = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dictolist_6tr.txt'
)
fmt_contig_list = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dictolist.txt'
)
# format get_contig_count_dict outputs
fmt_dict_contig_count = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/{fn_metat}-{fn_targets}-{fn_kallisto}-dict.json'
)
# format get_norm_factor outputs
fmt_sample_norm_factor = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/norm_factors/{fn_metat}-{fn_targets}-norm_factors'
    + '/{fn_metat}-{fn_targets}-{fn_kallisto}-norm_factor.txt'
)
# get_table_contig_norm_count
fmt_table_contig_norm_count = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/tables_norm_count/{fn_metat}-{fn_targets}-tables_norm_count'
    + '/{fn_metat}-{fn_targets}-{fn_kallisto}-table_norm_count.csv'
)
# get_table_KO_norm_count
fmt_table_KO_norm_count = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/tables_norm_count/{fn_metat}-{fn_targets}-tables_norm_count'
    + '/{fn_metat}-{fn_targets}-{fn_kallisto}-table_KO_norm_count.csv'
)
# combine_table_KO_norm_count
fmt_combine_table_KO_norm_count = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/tables_norm_count/combined'
    + '/{fn_targets}-table_samples_KOs_norm_count.csv'
)
# get_sum_counts
fmt_sum_counts = (
    config['snakemake_output_dir'] + '/count_sums/{fn_metat}-{fn_targets}-count_sums'
    + '/{fn_kallisto}-count_sum.txt'
)
# get_dict_contig_taxon
fmt_dict_contig_tax = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_tax/{fn_metat}-{fn_targets}-dict_contig_taxid.json'
)
# format taxon full table
fmt_df_contig_tax_full = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_tax/{fn_metat}-{fn_targets}-df_contig_full_taxonomy.csv'
)

fmt_sample_tidytable = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-{fn_kallisto}-tidy.csv'
)
fmt_sample_metadata = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-{fn_kallisto}-metadata.csv'
)

fmt_mergesample_tidytable = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-tidy.csv'
)
fmt_mergesample_metadata = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-metadata.csv'
)


fmt_mergesample_tidy_trim = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-tidy_trim.csv'
)
fmt_mergesample_dict_taxtrim_contigs = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-dict_taxtrim_contigs.json'
)
fmt_mergesample_tensor_taxtrim_tidy = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/{fn_metat}-{fn_targets}-tidys'
    + '/{fn_metat}-{fn_targets}-barnacle_tensor_tidy.csv'
)

fmt_mergeall_tidytable = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/merge_all'
    + '/{fn_targets}-tidy_all.csv'
)
fmt_mergeall_metadata = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/merge_all'
    + '/{fn_targets}-metadata.csv'
)

fmt_mergeall_dict_taxtrim_contigs = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/merge_all'
    + '/{fn_targets}-dict_taxtrim_contigs.csv'
)
fmt_mergeall_tensor_taxtrim_tidy = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/tidy_tables/merge_all'
    + '/{fn_targets}-barnacle_tensor_tidy.csv'
)

# fmt_mergeall_tidy_trim = (
#     config['snakemake_output_dir'] + '/{search_output_dir}'
#     + '/tidy_tables'
#     + '/{fn_targets}-tidy_all_trim.csv'
# )
# Write start files
fns_start = get_fns(fmt_start)
fn_other_start = config['snakemake_output_dir'] + '/start_files/start.txt'
fns_start.append(fn_other_start)
for fn in fns_start:
    if not os.path.exists(fn):
        # Make directories
        out_dir = os.path.split(fn)[0]
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # Write files
        with open(fn, 'w') as f:
            f.write('Started...')


# Scripts
sc_get_metat_dict = 'sc_get_metat_dict.py'
sc_get_column_sum = 'sc.get_column_sum.py'
# sc_get_metat_dict = '../repo-armbrust-metat-search/scripts/sc_get_metat_dict.py'


# =============================================================================
# Rule all output
# =============================================================================

dicts_done = get_fns(fmt_dict)
contig_list_done = get_fns(fmt_contig_list)
dicts_contig_count_done = get_fns_sample(fmt_dict_contig_count)
fns_sample_norm_factor = get_fns_sample(fmt_sample_norm_factor, col_check='sn_type_norm_factor')
fns_table_contig_norm_count = get_fns_sample(fmt_table_contig_norm_count, col_check='sn_type_norm_factor')
fns_table_KO_norm_count = get_fns_sample(fmt_table_KO_norm_count, col_check='sn_type_norm_factor')
fns_combine_table_KO_norm_count = get_fns_target(fmt_combine_table_KO_norm_count)
fns_sum_counts = get_fns_sample(fmt_sum_counts, col_check='sn_type_norm_factor')
fns_dict_contig_tax = get_fns(fmt_dict_contig_tax, col_check='fn_diamond')
fns_df_contig_tax_full = get_fns(fmt_df_contig_tax_full, col_check='fn_diamond')
fns_sample_tidytable = get_fns_sample(fmt_sample_tidytable, col_check='sn_type_norm_factor')
fns_sample_metadata = get_fns_sample(fmt_sample_metadata, col_check='sn_type_norm_factor')
fns_mergesample_tidytable = get_fns(fmt_mergesample_tidytable, col_check='fn_diamond')
fns_mergesample_tidy_trim = get_fns(fmt_mergesample_tidy_trim, col_check='fn_diamond')
fns_mergesample_tensor_taxtrim_tidy = get_fns(fmt_mergesample_tensor_taxtrim_tidy, col_check='fn_diamond')
fn_mergeall_tidytable = get_fns_target(fmt_mergeall_tidytable)
fn_mergeall_tidytable = get_fns_target(fmt_mergeall_tidytable)
fn_mergeall_tensor_taxtrim_tidy = get_fns_target(fmt_mergeall_tensor_taxtrim_tidy)

# =============================================================================
# Snake rules
# =============================================================================

rule all:
    input:
        fns_mergesample_tensor_taxtrim_tidy,
        fn_mergeall_tensor_taxtrim_tidy,

rule get_dict_KO_contig:
    input:
        fn_start = fmt_start
    output:
        fn_dict = fmt_dict
    params:
        script = sc_get_metat_dict,
        fn_metat_full = lambda w: get_fn_full(
            w.fn_metat, w.fn_targets, input_table,
            col_dir='path_metat', col_fn='fn_metat'
        ),
        fn_targets_full = lambda w: get_fn_full(
            w.fn_metat, w.fn_targets, input_table,
            col_dir='path_targets', col_fn='fn_targets'
        ),
        name_key = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'name_key'
        ),
        name_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'name_value'
        ),
        columns_line_number = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'columns_line_number'
        ),
        idx_key = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_key'
        ),
        idx_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_value'
        ),
    shell:
        """
        {params.script} \
            -m {params.fn_metat_full} \
            -t {params.fn_targets_full} \
            -k {params.name_key} \
            -vl {params.name_value}  \
            -c {params.columns_line_number} \
            -ik {params.idx_key} \
            -ivl {params.idx_value} \
            -o {output.fn_dict} \
            --verbose
        """
    

rule get_list_contigs:
    input:
        fn_dict = fmt_dict
    output:
        fn_contig_list = fmt_contig_list
    run:
        # Load dict
        with open(input.fn_dict, 'r') as f:
            dict_lists = json.load(f)
        
        # Merge and save
        with open(output.fn_contig_list, 'w') as fo:
            for _, names in dict_lists.items():
                for n_ in names:
                    n = re.split(r'_\d+$',n_)[0]  # Remove "_{digit}" aa reading frame from the contig name
                    fo.write(f'{n}\n')

rule get_list_contigs_6tr:
    input:
        fn_dict = fmt_dict
    output:
        fn_contig_list_6tr = fmt_contig_list_6tr,
    run:
        # Load dict
        with open(input.fn_dict, 'r') as f:
            dict_lists = json.load(f)
        
        # Merge and save
        with open(output.fn_contig_list_6tr, 'w') as fo6:
            for _, names in dict_lists.items():
                for n_ in names:
                    fo6.write(f'{n_}\n')


rule get_dict_contig_count:
    input:
        fn_contig_list = fmt_contig_list
    output:
        fn_dict_contig_count = fmt_dict_contig_count
    params:
        script = sc_get_metat_dict,
        fn_kallisto_full = lambda w: get_fn_kallisto_full(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        name_key = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'name_key_kallisto'
        ),
        name_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'name_value_kallisto'
        ),
        columns_line_number = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'col_line_kallisto'
        ),
        idx_key = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_key_kallisto'
        ),
        idx_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_value_kallisto'
        ),
    shell:
        """
        {params.script} \
            -m {params.fn_kallisto_full} \
            -mt {params.fn_kallisto_tar} \
            -t {input.fn_contig_list} \
            -k {params.name_key} \
            -vl {params.name_value}  \
            -c {params.columns_line_number} \
            -ik {params.idx_key} \
            -ivl {params.idx_value} \
            -o {output.fn_dict_contig_count} \
            --verbose
        """

rule get_norm_factor:
    input:
        fn_dict_contig_count = fmt_dict_contig_count
    output:
        fn_sample_norm_factor = fmt_sample_norm_factor
    params:
        fn_norm_factor = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'fn_norm_factor'
        ),
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        sn_type = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'sn_type_norm_factor'
        ),
        col_name = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'col_name_norm_factor'
        ),
    run:
        # Get kallisto counts file name
        fn_sample_counts = wildcards.fn_kallisto
        if params.fn_kallisto_tar:
            fn_sample_counts = params.fn_kallisto_tar  

        # Get sample name used in the norm factor table
        if params.sn_type == 'G1':
            sn = re.search(r'(?<=\.)S.+?um.*[A-Z]{1,2}(?=\.)', fn_sample_counts)[0]
            sn = re.sub('\.','',sn)
        elif params.sn_type == 'G2':
            # Get regex matches
            sns = re.findall(r'(?=(G.+?um.*[A-Z]{1,2}(?=\.)))', fn_sample_counts)
            # Pick the shortest match
            sns_len = [len(sn) for sn in sns]
            sns_len_min = min(sns_len)
            for i, l in enumerate(sns_len):
                if l == sns_len_min:
                    idx = i
            sn = sns[i]
        elif params.sn_type == 'D1':
            sn = re.search(r'S.+?(?=\.)', fn_sample_counts)[0]
        elif params.sn_type == 'G3PA':
            sn = re.search(r'G.+(?=\.unstranded)', fn_sample_counts)[0]
            sn += '.flash'
        elif params.sn_type == 'G3NS':
            sn = os.path.splitext(fn_sample_counts)[0]
        else:
            raise ValueError(
                f"""
                Sample name type not provided (sn_type column in file table)
                """
            )
        # elif params.sn_type == 'G2PA_DCM_RR':
        #     sn = re.search(r'^\w+(?=\.)', fn_sample_counts) + '.sam' 

        
        if sn:
            # Get the norm factor from the table
            df_norm = pd.read_csv(params.fn_norm_factor)
            try:
                norm_factor = df_norm.loc[
                    df_norm[params.col_name] == sn,
                    'NORM_FACTOR'
                ].values[0]
            except:
                raise ValueError(
                    f"""
                    The sample name matching failed. The filename to parse was: 
                    {fn_sample_counts}
                    The parsed sample name was:
                    {sn}
                    The dataframe of normalization factors was:
                    {df_norm}
                    The column to search for sample name was:
                    {params.col_name}
                    """
                ) 
        else:
            norm_factor = ''

        # Write norm factor to file
        with open(output.fn_sample_norm_factor, 'w') as f:
            f.write(str(norm_factor))


rule get_sum_counts:
    input:
        fn_start = fn_other_start
    output:
        fn_sum_counts = fmt_sum_counts
    params:
        script = sc_get_column_sum,
        fn_kallisto_full = lambda w: get_fn_kallisto_full(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        name_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'name_value_kallisto'
        ),
        columns_line_number = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'col_line_kallisto'
        ),
        idx_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_value_kallisto'
        ),
    shell:
        """
        {params.script} \
            -m {params.fn_kallisto_full} \
            -mt {params.fn_kallisto_tar} \
            -vl {params.name_value}  \
            -c {params.columns_line_number} \
            -ivl {params.idx_value} \
            -o {output.fn_sum_counts} \
            --verbose
        """        


rule get_table_contig_norm_count:
    input:
        fn_sample_norm_factor = fmt_sample_norm_factor,
        fn_dict_contig_count = fmt_dict_contig_count,
        fn_contig_list = fmt_contig_list,
    output:
        fn_table_norm_count = fmt_table_contig_norm_count
    run:
        # Load norm factor
        with open(input.fn_sample_norm_factor, 'r') as f:
            norm_factor = f.read()
        norm_factor = float(norm_factor)
        # Load dict
        with open(input.fn_dict_contig_count, 'r') as f:
            dict_contig_count = json.load(f)
        # Load list
        with open(input.fn_contig_list, 'r') as f:
            list_contigs = f.read().splitlines()

        # Write contigs with no kallisto alignment to separate file       
        bn, _ = os.path.splitext(output.fn_table_norm_count)
        fn_failed = bn + '-no_kallisto_align.txt'
        with open(output.fn_table_norm_count, 'w') as f, open(fn_failed, 'w') as ff:
            writer = csv.writer(f)
            # Write column names
            columns = ['contig', wildcards.fn_kallisto]
            writer.writerow(columns)
            # Write 
            for contig in list_contigs:
                # Normalize
                count = dict_contig_count.get(contig)
                if count:
                    count = float(count[0])
                    count *= norm_factor
                else:
                    # If not in the dictionary, note it in the failed file
                    count = 0
                    ff.write(f'{contig}\n')
                row = [contig, count]
                writer.writerow(row)

rule get_table_KO_norm_count:
    input:
        fn_dict_KO_contigs = fmt_dict,
        fn_dict_contig_count = fmt_dict_contig_count,
        fn_sample_norm_factor = fmt_sample_norm_factor,
        fn_sum_counts = fmt_sum_counts
    output:
        fn_table_KO_norm_count = fmt_table_KO_norm_count,
    params:
        fn_list_KOs = lambda w: get_fn_full(
            w.fn_metat, w.fn_targets, input_table,
            col_dir='path_targets', col_fn='fn_targets'
        ),
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        sn_type_parse_kallisto = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 
            'sn_type_parse_kallisto'
        ),
    run:
        # Load dict ko -> contigs
        with open(input.fn_dict_KO_contigs, 'r') as f:
            dict_KO_contigs = json.load(f)
        # Load dict contig -> count
        with open(input.fn_dict_contig_count, 'r') as f:
            dict_contig_count = json.load(f)
        # Load norm factor
        with open(input.fn_sample_norm_factor, 'r') as f:
            norm_factor = f.read()
        norm_factor = float(norm_factor)
        # Load sum counts
        with open(input.fn_sum_counts, 'r') as f:
            sum_counts = f.read()
        sum_counts = float(sum_counts)
        # Convert to dict ko -> norm_counts, summing across all contigs in the KO
        dict_KO_norm_count = defaultdict(lambda: 0)
        for KO, contigs in dict_KO_contigs.items():
            for c_ in contigs:
                c = re.split(r'_\d+$',c_)[0]  # Remove "_{digit}" aa reading frame from the contig name
                count = dict_contig_count.get(c)
                if count:
                    count = float(count[0])
                    count *= norm_factor
                else:
                    # If not in the dictionary
                    count = 0
                dict_KO_norm_count[KO] += count
        # Load KO list
        with open(params.fn_list_KOs) as f:
            list_KOs = f.read().splitlines()
        # Get the kallisto filename
        fn_sample_counts = wildcards.fn_kallisto
        if params.fn_kallisto_tar:
            fn_sample_counts = params.fn_kallisto_tar  
        # Start the row
        columns = ['fn_KO', 'fn_sample_counts']
        row = [wildcards.fn_metat, fn_sample_counts]
        # Parse the sample name info
        columns += parse_fn_kallisto_sn(
            fn_sample_counts,
            get_columns=True
        )
        row += parse_fn_kallisto_sn(
            fn_sample_counts, 
            params.sn_type_parse_kallisto
        )
        # Add sum counts for the sample
        columns += ['sample_sum_counts']
        row += [sum_counts * norm_factor]
        # Build row in order by ko list 
        columns += list_KOs
        row += [dict_KO_norm_count[KO] for KO in list_KOs]
        # Write to file
        with open(output.fn_table_KO_norm_count, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            writer.writerow(row)


rule combine_table_KO_norm_count:
    input:
        fns_table_KO_norm_count = fns_table_KO_norm_count,
    output:
        fn_combine_table_KO_norm_count = fmt_combine_table_KO_norm_count
    run:
        with open(output.fn_combine_table_KO_norm_count, 'w') as fw:
            i = 0
            for fn in input.fns_table_KO_norm_count:
                with open(fn, 'r') as fr:
                    columns = fr.readline()
                    if i == 0:
                        fw.write(columns)
                    row = fr.readline()
                    fw.write(row)
                i += 1




rule get_dict_contig_taxon:
    input:
        fn_contig_list_6tr = fmt_contig_list_6tr,
        fn_contig_list = fmt_contig_list
    output:
        fn_dict_contig_tax = fmt_dict_contig_tax
    params:
        script = sc_get_metat_dict,
        fn_diamond = lambda w: get_fn_full(
            w.fn_metat, w.fn_targets, input_table,
            col_dir='dir_diamond', col_fn='fn_diamond'
        ),
        idx_key = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_key_diamond'
        ),
        idx_value = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 'idx_value_diamond'
        ),
        fn_list_contigs = lambda w, input: get_fn_list_contigs_diamond(
            w.fn_metat, w.fn_targets, input_table,
            fn_6tr=input.fn_contig_list_6tr, 
            fn_std=input.fn_contig_list
        ),
    shell:
        """
        {params.script} \
            -m {params.fn_diamond} \
            -t {params.fn_list_contigs} \
            -ik {params.idx_key} \
            -ivl {params.idx_value} \
            -o {output.fn_dict_contig_tax} \
            --verbose
        """

rule get_table_full_lineage:
    input:
        fn_dict_contig_tax = fmt_dict_contig_tax,
        fn_contig_list_6tr = fmt_contig_list_6tr,
        fn_contig_list = fmt_contig_list,
    output:
        fn_df_contig_tax_full = fmt_df_contig_tax_full
    params:
        fn_list_contigs = lambda w, input: get_fn_list_contigs_diamond(
            w.fn_metat, w.fn_targets, input_table,
            fn_6tr=input.fn_contig_list_6tr, 
            fn_std=input.fn_contig_list
        )
    run:
        # Load list
        with open(params.fn_list_contigs, 'r') as f:
            list_contigs = f.read().splitlines()
            
        # Load dict
        with open(input.fn_dict_contig_tax, 'r') as f:
            dict_contig_taxid = json.load(f)

        # Get database        
        ncbi = NCBITaxa()
        # TODO build in a function to check if database needs to be updated

        bn, _ = os.path.splitext(output.fn_df_contig_tax_full)
        fn_failed = bn + '-failed_to_fetch.txt'
        with open(output.fn_df_contig_tax_full, 'w') as f, open(fn_failed, 'w') as ff:
            writer = csv.writer(f)
            # Write columns
            columns = ['contig']
            for rank in config['desired_ranks']:
                columns.append(rank + '_name')
                columns.append(rank + '_taxid')
            writer.writerow(columns)
            # Write a row of ranks for each contig
            dict_rank_taxid = {}
            for contig in list_contigs:
                # Get lineage for contig
                taxid = dict_contig_taxid.get(contig)
                if taxid:
                    taxid = taxid[0]
                    if not int(taxid) == 0:
                        try:
                            lineage = ncbi.get_lineage(taxid)
                            # Get ranks for lineage
                            dict_taxid_rank = ncbi.get_rank(lineage)
                            dict_rank_taxid = {v:k for k, v in dict_taxid_rank.items()}
                            # Get names for lineage
                            names = ncbi.translate_to_names(lineage)
                            dict_taxid_name = dict(zip(lineage, names))
                        except:
                            ff.write(f'{taxid}\n')
                            pass
                # Get a list of names and taxids for the lineage ordered by rank
                tax_full = [contig]
                for rank in config['desired_ranks']:
                    taxid = dict_rank_taxid.get(rank)
                    if taxid:
                        name = dict_taxid_name[taxid]
                    else:
                        name = None
                    tax_full.append(name)
                    tax_full.append(taxid)
                # Write to file
                writer.writerow(tax_full)
            
rule get_sample_tidytable:
    input:
        fn_dict_KO_contigs = fmt_dict,
        fn_dict_contig_tax = fmt_dict_contig_tax,
        fn_dict_contig_count = fmt_dict_contig_count
    output:
        fn_sample_tidytable = fmt_sample_tidytable,
        # fn_sample_metadata = fmt_sample_metadata
    params:
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        sn_type_parse_kallisto = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 
            'sn_type_parse_kallisto'
        ),
        type_contig_tax = lambda w: get_input_table_value(
                w.fn_metat, w.fn_targets, input_table, 'type_diamond_contig_name'
            )
    run:
        # Load kos
        with open(input.fn_dict_KO_contigs, 'r') as f:
            dict_ko_contigs = json.load(f)
        # Invert dict
        dict_contig_ko = {}
        for ko, contigs in dict_ko_contigs.items():
            for c in contigs:
                dict_contig_ko[c] = ko
        # Load taxa
        with open(input.fn_dict_contig_tax, 'r') as f:
            dict_contig_tax = json.load(f)
        # Load counts
        with open(input.fn_dict_contig_count, 'r') as f:
            dict_contig_count = json.load(f)
        # Get the kallisto filename
        fn_sample_counts = wildcards.fn_kallisto
        if params.fn_kallisto_tar:
            fn_sample_counts = params.fn_kallisto_tar 
        # # Get sample metdata columns
        # sample_info_columns = parse_fn_kallisto_sn('', get_columns=True)
        # # Get sample metadata
        # sample_info = parse_fn_kallisto_sn(
        #     fn_sample_counts, 
        #     params.sn_type_parse_kallisto
        # )
        # # Get sample name without replicate
        # namenorep = ''
        # for col, info in zip(sample_info_columns, sample_info):
        #     if col != 'rep':
        #         namenorep += f"{str(row[c])}-"
        # assm_sample = namenorep[:-1]
        # # Add info to rows
        # sample_info_columns += ['fn_sample_counts', 'assm_sample']
        # sample_info += [fn_sample_counts, assm_sample]
        # # Write to metadata file
        # with open(output.fn_sample_metadata, 'w') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(sample_info_columns)
        #     writer.writerow(sample_info)        
        # Get sample plus assembly identifier
        # dict_si = dict(zip(sample_info_columns, sample_info))
        # assm_sample = dict_si['assembly'] + '-' + dict_si['sample']
        # add identifiers to columns and metadata 
        # sample_info_columns += ['fn_sample_counts', 'assm_sample']
        # sample_info += [fn_sample_counts, assm_sample]
        # sample_info_columns.append('fn_sample_counts')
        # sample_info.append(fn_sample_counts)
        # # Write to metadata file
        # with open(output.fn_sample_metadata, 'w') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(sample_info_columns)
        #     writer.writerow(sample_info)
        # Get columns list
        columns = ['contig','fn_sample_counts', 'KO','taxon','estcounts']
        rep = dict(zip(sample_info_columns, sample_info))['rep']
        # Write lines
        with open(output.fn_sample_tidytable, 'w') as fo:
            writer = csv.writer(fo)
            writer.writerow(columns)
            countsum = 0
            n_tax = 0
            for c, ko in dict_contig_ko.items():
                c_ = re.sub(r'_\d+$','',c)      
                ctax = c if params.type_contig_tax else c_
                tax = dict_contig_tax.get(ctax)
                tax = tax[0] if tax else 0
                if tax:
                    n_tax += 1
                # tax = dict_contig_tax[c].get(c)
                # tax = tax[0] if tax else 0
                count = dict_contig_count.get(c_)
                count = count[0] if count else 0
                countsum += float(count)
                row = [c, fn_sample_counts, ko, tax, count]
                writer.writerow(row)
        if countsum == 0:
            raise ValueError(
                    f"""
                    Counts were not read from {input.fn_dict_contig_count}, or there are no counts.
                    """
                ) 
        if n_tax == 0:
            raise ValueError(
                    f"""
                    Taxa were not read from {input.fn_dict_contig_count}, or there are no taxonomic annotations.
                    """
                ) 

rule get_sample_metadata:
    input:
        fn_dict_KO_contigs = fmt_dict,
        fn_dict_contig_tax = fmt_dict_contig_tax,
        fn_dict_contig_count = fmt_dict_contig_count
    output:
        fn_sample_metadata = fmt_sample_metadata
    params:
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        sn_type_parse_kallisto = lambda w: get_input_table_value(
            w.fn_metat, w.fn_targets, input_table, 
            'sn_type_parse_kallisto'
        ),
        type_contig_tax = lambda w: get_input_table_value(
                w.fn_metat, w.fn_targets, input_table, 'type_diamond_contig_name'
            )
    run:
        # Get the kallisto filename
        fn_sample_counts = wildcards.fn_kallisto
        if params.fn_kallisto_tar:
            fn_sample_counts = params.fn_kallisto_tar 
        # Get sample metdata columns
        sample_info_columns = parse_fn_kallisto_sn('', get_columns=True)
        # Get sample metadata
        sample_info = parse_fn_kallisto_sn(
            fn_sample_counts, 
            params.sn_type_parse_kallisto
        )
        # Get sample name without replicate
        namenorep = ''
        for col, info in zip(sample_info_columns, sample_info):
            if col != 'rep':
                if info:
                    namenorep += f"{info}-"
        assm_sample = namenorep[:-1]
        # Add info to rows
        sample_info_columns += ['fn_sample_counts', 'assm_sample']
        sample_info += [fn_sample_counts, assm_sample]
        # Write to metadata file
        with open(output.fn_sample_metadata, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(sample_info_columns)
            writer.writerow(sample_info)




rule merge_sample_tidy_tables:
    input:
        fns_sample_tidytable = lambda w: get_fns_row_sample(
            fmt_sample_tidytable,
            w.fn_metat, w.fn_targets, input_table
        ),
        # fn_sample_metadata = lambda w: get_fns_row_sample(
        #     fmt_sample_metadata,
        #     w.fn_metat, w.fn_targets, input_table
        # ),
    output:
        fn_mergesample_tidytable = fmt_mergesample_tidytable,
        fn_mergesample_metadata = fmt_mergesample_metadata
    run:
        # Write table
        with open(output.fn_mergesample_tidytable, 'w') as fw:
            i = 0
            for fn in input.fns_sample_tidytable:
                with open(fn, 'r') as fr:
                    columns = fr.readline()
                    if i == 0:
                        fw.write(columns)
                    for row in fr:
                        fw.write(row)
                i += 1   
        # # Write metadata     
        # with open(output.fn_mergesample_metadata, 'w') as fw:
        #     i = 0
        #     for fn in input.fn_sample_metadata:
        #         with open(fn, 'r') as fr:
        #             columns = fr.readline()
        #             row = fr.readline()
        #             if i == 0:
        #                 fw.write(columns)
        #             fw.write(row)
        #         i += 1

rule merge_all_metadata:
    input:
        fns_sample_metadata = fns_sample_metadata,
    output:
        fn_mergeall_metadata = fmt_mergeall_metadata,
    run:
        # Write metadata     
        with open(output.fn_mergeall_metadata, 'w') as fw:
            i = 0
            for fn in input.fns_sample_metadata:
                with open(fn, 'r') as fr:
                    columns = fr.readline()
                    row = fr.readline()
                    if i == 0:
                        fw.write(columns)
                    fw.write(row)
                i += 1 

rule merge_all_tidy_tables:
    input:
        fns_sample_tidytable = fns_sample_tidytable,
        # fns_sample_metadata = fns_sample_metadata,
    output:
        fn_mergeall_tidytable = fmt_mergeall_tidytable,
        # fn_mergeall_metadata = fmt_mergeall_metadata,
    run:
        # Write table
        with open(output.fn_mergeall_tidytable, 'w') as fw:
            i = 0
            for fn in input.fns_sample_tidytable:
                with open(fn, 'r') as fr:
                    columns = fr.readline()
                    if i == 0:
                        fw.write(columns)
                    for row in fr:
                        fw.write(row)
                i += 1   
        # # Write metadata     
        # with open(output.fn_mergeall_metadata, 'w') as fw:
        #     i = 0
        #     for fn in input.fns_sample_metadata:
        #         with open(fn, 'r') as fr:
        #             columns = fr.readline()
        #             row = fr.readline()
        #             if i == 0:
        #                 fw.write(columns)
        #             fw.write(row)
        #         i += 1 


rule trim_sample_taxon_trees_build_tidy_tensor:
    input:
        fn_mergesample_tidytable = fmt_mergesample_tidytable,
        fn_mergeall_metadata = fmt_mergeall_metadata,
        # fn_mergesample_metadata = fmt_mergesample_metadata
    output:
        fn_mergesample_dict_taxtrim_contigs = fmt_mergesample_dict_taxtrim_contigs,
        fn_mergesample_tensor_taxtrim_tidy = fmt_mergesample_tensor_taxtrim_tidy
        # fn_mergesample_tidy_trim = fmt_mergesample_tidy_trim
    params:
        filtfunc = config['tree_trim']['filtfunc'],
        thresh = config['tree_trim']['thresh'],
        minsamples = config['tree_trim']['minsamples'],
    run:
        # Build tree
        t = TreeTrim(
            input.fn_mergesample_tidytable,
            col_contig='contig', 
            col_ko='KO',
            col_taxon='taxon',
            col_estcounts='estcounts',
            col_sample='fn_sample_counts'
        )
        # Trim tree
        t.trim_tree(
            filt_func_name=params.filtfunc,
            thresh=int(params.thresh), 
            minsamples=int(params.minsamples)
        )
        # save dict taxtrim contig
        with open(output.fn_mergesample_dict_taxtrim_contigs, 'w') as f:
            json.dump(
                t.dict_taxtrim_contigs, 
                f, 
                sort_keys=True, 
                indent=4, 
                separators=(',', ': ')
            )
        # Get dict sample rep
        metadata = pd.read_csv(input.fn_mergeall_metadata)
        dict_sam_namenorep = dict(zip(metadata['fn_sample_counts'], metadata['assm_sample']))
        dict_sam_rep = dict(zip(metadata['fn_sample_counts'], metadata['rep']))
        # cols_namenorep = metadata.columns[
        #     (metadata.columns != 'fn_sample_counts')
        #     & (metadata.columns != 'rep')
        # ]
        # dict_sam_rep = {}
        # dict_sam_namenorep = {}
        # for _, row in metadata.iterrows():
        #     dict_sam_rep[row['fn_sample_counts']] = row['rep']
        #     namenorep = ''
        #     for c in cols_namenorep:
        #         namenorep += f"{str(row[c])}-"
        #     dict_sam_namenorep[row['fn_sample_counts']] = namenorep[:-1]
        # dict_sam_rep = dict(zip(
        #     metadata['assm_sample'], metadata['rep']
        # ))
        # Build tensor sample x gene x taxon as tidytable
        # dict_taxtrim_ko_sample = defaultdict(lambda: defaultdict(dict))
        columns = ['assm_sample','KO','taxon_trim','estcounts','rep']
        with open(output.fn_mergesample_tensor_taxtrim_tidy, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            for taxtrim, d1 in t.dict_taxtrim_ko_sam_estcounts.items():
                for ko, d2 in d1.items():
                    for sam, ec in d2.items():
                        rep = dict_sam_rep[sam]
                        namenorep = dict_sam_namenorep[sam]
                        row = [namenorep, ko, taxtrim, ec, rep]
                        writer.writerow(row)

        # # Add trimmed taxa column to tidy table
        # with open(input.fn_mergesample_tidytable, 'r') as fr:
        #     with open(output.fn_mergesample_tidy_trim, 'w') as fw:
        #         reader = csv.DictReader(fr)
        #         columns = reader.fieldnames + ['taxon_trim']
        #         writer = csv.DictWriter(fw, delimiter=',', fieldnames=columns)
        #         # Write columns line
        #         writer.writeheader()
        #         # write rest of lines
        #         for row in reader:
        #             contig = row['contig']
        #             taxtrim = t.dict_contig_taxtrim.get(contig)
        #             taxtrim = taxtrim if taxtrim else 0
        #             row['taxon_trim'] = taxtrim
        #             writer.writerow(row)


rule trim_merge_all_taxon_tree_build_tidy_tensor:
    input:
        fn_mergeall_tidytable = fmt_mergeall_tidytable,
        fn_mergeall_metadata = fmt_mergeall_metadata,
    output:
        fn_mergeall_dict_taxtrim_contigs = fmt_mergeall_dict_taxtrim_contigs,
        fn_mergeall_tensor_taxtrim_tidy = fmt_mergeall_tensor_taxtrim_tidy,
    params:
        filtfunc = config['tree_trim_merge']['filtfunc'],
        thresh = config['tree_trim_merge']['thresh'],
        minsamples = config['tree_trim_merge']['minsamples'],
        minbatches = config['tree_trim_merge']['minbatches'],
    run:
        # Get sample batch dict
        metadata = pd.read_csv(input.fn_mergeall_metadata)
        dict_sample_batch = dict(zip(
            metadata['fn_sample_counts'],
            metadata['assembly']
        ))

        # Build tree
        t = TreeTrim(
            input.fn_mergeall_tidytable,
            col_contig='contig', 
            col_ko='KO',
            col_taxon='taxon',
            col_estcounts='estcounts',
            col_sample='fn_sample_counts'
        )
        # Trim tree
        t.trim_tree(
            filt_func_name=params.filtfunc,
            thresh=int(params.thresh),
            dict_sample_batch=dict_sample_batch,
            minsamples=int(params.minsamples),
            minbatches=int(params.minbatches),
        )
        # save dict taxtrim contig
        with open(output.fn_mergeall_dict_taxtrim_contigs, 'w') as f:
            json.dump(
                t.dict_taxtrim_contigs, 
                f, 
                sort_keys=True, 
                indent=4, 
                separators=(',', ': ')
            )
        # Get dict sample rep
        dict_sam_namenorep = dict(zip(metadata['fn_sample_counts'], metadata['assm_sample']))
        dict_sam_rep = dict(zip(metadata['fn_sample_counts'], metadata['rep']))
        # cols_namenorep = metadata.columns[
        #     (metadata.columns != 'fn_sample_counts')
        #     & (metadata.columns != 'rep')
        # ]
        # dict_sam_rep = {}
        # dict_sam_namenorep = {}
        # for _, row in metadata.iterrows():
        #     dict_sam_rep[row['fn_sample_counts']] = row['rep']
        #     namenorep = ''
        #     for c in cols_namenorep:
        #         namenorep += f"{str(row[c])}-"
        #     dict_sam_namenorep[row['fn_sample_counts']] = namenorep[:-1]
        # dict_sam_rep = dict(zip(
        #     metadata['assm_sample'], metadata['rep']
        # ))
        # Build tensor sample x gene x taxon as tidytable
        # dict_taxtrim_ko_sample = defaultdict(lambda: defaultdict(dict))
        columns = ['assm_sample','KO','taxon_trim','estcounts','rep']
        with open(output.fn_mergeall_tensor_taxtrim_tidy, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            for taxtrim, d1 in t.dict_taxtrim_ko_sam_estcounts.items():
                for ko, d2 in d1.items():
                    for sam, ec in d2.items():
                        rep = dict_sam_rep[sam]
                        namenorep = dict_sam_namenorep[sam]
                        row = [namenorep, ko, taxtrim, ec, rep]
                        writer.writerow(row)

