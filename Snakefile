#!/usr/bin/env -S snakemake -s

"""
Pipeline to extract data from the Armbrust Gradients metaT dataset

"""

# =============================================================================
# Imports
# =============================================================================

import os
import json
import math
import glob
import tarfile
import fnmatch
import pandas as pd

# =============================================================================
# Functions
# =============================================================================

def get_tar_names(fn_tar, re_tar, search_output_dir):
    # Make tar name dir
    out_dir = config['snakemake_output_dir'] + f'/{search_output_dir}/tar_names'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Make filename
    fn_tar_names = os.path.split(fn_tar)[1]
    fn_tar_names = re.sub('.tar.gz','-tar_names.txt', fn_tar_names)
    fn_tar_names_full = f'{out_dir}/{fn_tar_names}'
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

def clip_path(fn):
    return os.path.split(fn)[1]

def get_fns(fmt):
    fns = []
    for index, row in input_table.iterrows():
        # Build and append filename
        fns.append(fmt.format(
            search_output_dir=row.search_output_dir,
            fn_metat=row.fn_metat, 
            fn_targets=row.fn_targets,
        ))
    return fns

def get_fns_sample(fmt):
    fns = []
    for index, row in input_table.iterrows():
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
                    row.search_output_dir
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

def get_fn_metat_full(fn_metat, fn_targets, input_table):
    path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_metat')
    if path:
        return path + f'/{fn_metat}'
    else:
        return fn_metat

def get_fn_targets_full(fn_metat, fn_targets, input_table):
    path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_targets')
    if path:
        return path + f'/{fn_targets}'
    else:
        return fn_targets

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


# =============================================================================
# Setup
# =============================================================================


input_table = pd.read_csv(
    config['input_table'],
    keep_default_na=False
)
# # Replace nan with empty string
# input_table.fillna('', inplace=True) 
# int_columns = [
#     'columns_line_number',
#     'idx_key',
#     'idx_value',
# ]
# for col in int_columns:
#     input_table[col] = input_table[col].astype(int)


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
fmt_contig_list = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dictolist.txt'
)

# format get_contig_count_dict outputs
fmt_dict_contig_count = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_count/{fn_metat}-{fn_targets}-{fn_kallisto}-dict.json'
)

# Write start files
fns_start = get_fns(fmt_start)
for fn in fns_start:
    # Make directories
    out_dir = os.path.split(fn)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Write files
    with open(fn, 'w') as f:
        f.write('Started...')


# =============================================================================
# Rule all output
# =============================================================================

dicts_done = get_fns(fmt_dict)
contig_list_done = get_fns(fmt_contig_list)
dicts_contig_count_done = get_fns_sample(fmt_dict_contig_count)
# =============================================================================
# Snake rules
# =============================================================================

rule all:
    input:
        dicts_contig_count_done


rule get_KO_contig_dict:
    input:
        fn_start = fmt_start
    output:
        fn_dict = fmt_dict
    params:
        fn_metat_full = lambda w: get_fn_metat_full(
            w.fn_metat, w.fn_targets, input_table
        ),
        fn_targets_full = lambda w: get_fn_targets_full(
            w.fn_metat, w.fn_targets, input_table
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
        sc_get_metat_dict.py \
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
    

rule get_contig_list:
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


rule get_contig_count_dict:
    input:
        fn_contig_list = fmt_contig_list
    output:
        fn_dict_contig_count = fmt_dict_contig_count
    params:
        fn_kallisto_full = lambda w: get_fn_kallisto_full(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        fn_kallisto_tar = lambda w: get_fn_kallisto_tar(
            w.fn_metat, w.fn_targets, w.fn_kallisto, input_table
        ),
        fn_targets_full = lambda w: get_fn_targets_full(
            w.fn_metat, w.fn_targets, input_table
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
        sc_get_metat_dict.py \
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
    