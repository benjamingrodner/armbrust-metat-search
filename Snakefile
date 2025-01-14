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
os.environ['XDG_DATA_HOME'] = config['XDG_DATA_HOME'] or os.environ['PWD']  # Dir to write the ncbitaxa database (written into 'ete' added to path)
from ete4 import NCBITaxa


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
# format get_norm_factor outputs
fmt_dict_contig_tax = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_tax/{fn_metat}-{fn_targets}-dict_contig_taxid.json'
)
# format taxon full table
fmt_df_contig_tax_full = (
    config['snakemake_output_dir'] + '/{search_output_dir}'
    + '/dicts_contig_tax/{fn_metat}-{fn_targets}-df_contig_full_taxonomy.csv'
)

# Write start files
fns_start = get_fns(fmt_start)
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
# sc_get_metat_dict = '../repo-armbrust-metat-search/scripts/sc_get_metat_dict.py'


# =============================================================================
# Rule all output
# =============================================================================

dicts_done = get_fns(fmt_dict)
contig_list_done = get_fns(fmt_contig_list)
dicts_contig_count_done = get_fns_sample(fmt_dict_contig_count)
fns_sample_norm_factor = get_fns_sample(fmt_sample_norm_factor)
fns_dict_contig_tax = get_fns(fmt_dict_contig_tax, col_check='fn_diamond')
fns_df_contig_tax_full = get_fns(fmt_df_contig_tax_full, col_check='fn_diamond')

# =============================================================================
# Snake rules
# =============================================================================

rule all:
    input:
        fns_sample_norm_factor,
        fns_df_contig_tax_full


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
            sn = ''
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
                print('df_norm')
                print(df_norm)
                print('\ncol_name')
                print(params.col_name)
                print('\nsn')
                print(sn)
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
            
                

