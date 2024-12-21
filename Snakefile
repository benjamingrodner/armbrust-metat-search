#!/usr/bin/env -S snakemake -s

"""
Pipeline to extract data from the Armbrust Gradients metaT dataset

"""

# =============================================================================
# Imports
# =============================================================================

import pandas as pd

# =============================================================================
# Functions
# =============================================================================

def get_input_table():
    input_table = pd.read_csv(config['input_table'])
    # input_table.columns = config['input_table_cols']
    return input_table
    
def get_fns(fmt):
    fns = []
    for index, row in input_table.iterrows():
        # # Get file basenames
        # bns = []
        # for fn_full in [row.fn_metat, args.fn_targets]:
        #     fn = os.path.split(fn_full)[1]
        #     bns.append(os.path.splitext(fn)[0])
        # Build and append filename
        fns.append(fmt.format(
            search_output_dir=row.search_output_dir,
            fn_metat=row.fn_metat, 
            fn_targets=row.fn_targets,
            name_key=row.name_key,
            name_value=row.name_value,
        ))
    return fns

def get_columns_line_number(fn_metat, input_table):
    return input_table.loc[
        input_table.fn_metat == fn_metat, 
        'columns_line_number'
    ].values[0]


# =============================================================================
# Setup
# =============================================================================


input_table = get_input_table()


# format start files
start_fmt = (
    config['snakemake_output_dir'] + '/start_files/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-{name_key}-{name_value}.txt'
)

# format get_metat_dict 
dict_fmt = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/{fn_metat}-{fn_targets}-dict-{name_key}-{name_value}.json'
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
        columns_line_number = lambda w: get_columns_line_number(w.fn_metat, input_table)
    shell:
        """
        sc_get_metat_dict.py \
            -m {wildcards.fn_metat} \
            -t {wildcards.fn_targets} \
            -k {wildcards.name_key} \
            -vl {wildcards.name_value}  \
            -c {params.columns_line_number} \
            -o {output.dict_fn} \
            --verbose
        """
    
