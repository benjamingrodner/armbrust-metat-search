#!/usr/bin/env -S snakemake -s

"""
Pipeline to extract data from the Armbrust Gradients metaT dataset

"""

# =============================================================================
# Imports
# =============================================================================

import math
import pandas as pd

# =============================================================================
# Functions
# =============================================================================
    
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

def get_input_table_value(fn_metat, fn_targets, input_table, column):
    return input_table.loc[
        (
            (input_table.fn_metat == fn_metat) 
            & (input_table.fn_targets == fn_targets)
        ), 
        column
    ].values[0]

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


# =============================================================================
# Setup
# =============================================================================


input_table = pd.read_csv(
    config['input_table'],
    keep_default_na=False
)
# # Replace nan with empty string
# input_table.fillna(0, inplace=True) 
# int_columns = [
#     'columns_line_number',
#     'idx_key',
#     'idx_value',
# ]
# for col in int_columns:
#     input_table[col] = input_table[col].astype(int)


# format start files
start_fmt = (
    config['snakemake_output_dir'] + '/start_files/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}.txt'
)

# format get_metat_dict 
dict_fmt = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dict.json'
)

# Write start files
start_fns = get_fns(start_fmt)
for fn in start_fns:
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

dicts_done = get_fns(dict_fmt)

# =============================================================================
# Snake rules
# =============================================================================

rule all:
    input:
        dicts_done


rule get_metat_dict:
    input:
        start_file = start_fmt
    output:
        dict_fn = dict_fmt
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
            -o {output.dict_fn} \
            --verbose
        """
    
