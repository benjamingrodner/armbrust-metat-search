#!/usr/bin/env python

"""
Script Name: sc_get_ko_contig_dicts.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2024-12-19
Description:
    Given a file with a list of Kegg orthologies
    Extract contigs from the gradients dataset with those kegg orthologies

Usage:
    python sc_get_ko_contig_dicts.py -m <fn_metat> -t <fn_targets> -k <name_key> -v <name_value> -c <columns_line_number> -ik <idx_key> -ivl <idx_value> -o <output_fn> --verbose

"""

###################################################################################################
# Imports
###################################################################################################

from collections import defaultdict
from time import time
import argparse
import tarfile
import gzip
import json
import sys
import os
import re
import io

###################################################################################################
# Functions
###################################################################################################
def save_dict(d, fn):
    """
    Write a human readable json file

    Args:
        d (dict): Dictionary to write.
        fn (str): Json output filename.
    """
    # Make folder if not existing
    od = os.path.split(fn)[0]
    if not os.path.exists(od):
        os.makedirs(od)
        print(f"Made dir: {od}")
    with open(fn, 'w') as f:
        json.dump(
            d, 
            f, 
            sort_keys=True, 
            indent=4, 
            separators=(',', ': ')
        )
    print(f"Wrote: {fn}")


def get_output_fn(args):
    """
    Build the output filename.

    Args:
        args (argparse.ArgumentParser): Command line arguments.

    Returns:
        str: Output filename
    """
    # bns = []
    # for fn_full in [args.fn_metat, args.fn_targets]:
    #     fn = os.path.split(fn_full)[1]
    #     bns.append(os.path.splitext(fn)[0])
    return(
        'output/' 
        + args.fn_metat + '-' + args.fn_targets 
        + '.dict-' + args.name_key + '-' + args.name_value 
        + '.json'
    )


def split_line(line, fn, fn_tar=''):
    """
    Split line string read from file

    Args:
        line (str): Unparsed line from the file.
        fn (str): Path to file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
    Returns:
        list: Split line.
    """
    fspl, ext = os.path.splitext(fn)
    # Check if gzipped
    if ext == '.gz':
        ext = os.path.splitext(fspl)[1]
    # Adjust if it's a tarball
    if fn_tar:
        ext = os.path.splitext(fn_tar)[1]
        line = line.decode('utf-8')
    # Separate based on file type
    if ext == '.csv':
        return re.split(r',', line)
    if ext == '.tsv':
        return re.split(r'\t', line)
    if ext == '.tab':
        return re.split(r'\s+', line)
    

def add_to_dict(line, fn, fn_tar, set_search, idx_key, idx_value, dict_metat):
    """
    Extract key and value from a line and add to the dictionary

    Args:
        line (str): Unparsed line from the file.
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
        set_search (set): Set of keys to add to the dictionary.
        idx_key (int): Index of the key in the parsed line (parsed by ',' or '\s').
        idx_value (int): Index of the mapped value in the parsed line (parsed by ',' or '\s').
        dict_metat (defaultdict(list)): Existing dictionary.

    Returns:
        defaultdict(list): Updated dictionary.
    """
    # Separate the line by spaces (\s) or commas
    l = split_line(line, fn, fn_tar)
    if len(l) > max(idx_key, idx_value):
        # Get key and remove extra quotes
        key = l[idx_key].replace('"','').replace("'", "")
        # Search target keys
        if key in set_search:
            # Add value to dict
            value = l[idx_value].replace('"','').replace("'", "")
            dict_metat[key].append(value)
    return dict_metat


def open_file(fn, fn_tar=''):
    """
    Open text file if gzipped, or not

    Args:
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.

    Returns:
        _io.TextIOWrapper: Metat file.
    """
    # Get file type
    ext = os.path.splitext(fn)[1]
    # Unzip to read
    if not fn_tar:
        if ext == '.gz':
            return gzip.open(fn, 'rt')
        else:
            return open(fn, 'r')
    else:
        tarf = tarfile.open(fn, "r|gz")
        for t in tarf:
            if fn_tar in t.name:
                return tarf.extractfile(t)


def get_metat_dict(fn, fn_tar, set_search, idx_key, idx_value):
    """
    Collect all contigs with given Kegg Orthologies
    SET UP to parse armbrust-metat folder files as of 2024-12-19

    Args:
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
        set_search (set): Set of keys to add to the dictionary.
        idx_key (int): Index of the key in the parsed line (parsed by ',' or '\s').
        idx_value (int): Index of the mapped value in the parsed line (parsed by ',' or '\s').

    Returns:
        defaultdict(list): Map of metat values key -> values.
    """
    dict_metat = defaultdict(list)
    f = open_file(fn, fn_tar)
    for line in f:
        dict_metat = add_to_dict(
            line, fn, fn_tar,
            set_search, idx_key, idx_value,
            dict_metat
        )
    f.close()
    return dict_metat


def get_key_value_idx(fn, fn_tar, name_key, name_value, cln=1):
    """
    Get indices for keys and values in parsed metat file.

    Args:
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
        names_key (str): Column name in metat file to be used for dictionary keys.
        names_value (str): Column name in metat file to be used for dictionary values.

    Returns:
        int: Index of the key in the parsed line (parsed by ',' or '\s').
        int: Index of the mapped value in the parsed line (parsed by ',' or '\s').
    """
    f = open_file(fn, fn_tar)
    for _ in range(cln):
        line = f.readline()
    f.close()
    l = split_line(line, fn, fn_tar)
    for idx, column in enumerate(l):
        # Remove extra quotes
        column = column.replace('"','').replace("'", "")
        if column == name_key:
            idx_key = idx
        elif column == name_value:
            idx_value = idx
    try:
        idx_key
    except NameError:
        raise ValueError(
            f"The provided argument -k {name_key} did not map to column names in the file.\n\
            Here is the split line (line number {cln} in the file) used to search for column names:\n\
            {l}\n\
            To change the column line, specify -c <columns_line_number> in the arguments (counting starts from 1 not 0)"
        )
    try:
        idx_value
    except NameError:
        raise ValueError(
            f"The provided argument -vl {name_value} did not map to column names in the file.\n\
            Here is the split line (line number {cln} in the file) used to search for column names:\n\
            {l}\n\
            To change the column line, specify -c <columns_line_number> in the arguments (counting starts from 1 not 0)"
        )
    return(idx_key, idx_value)

    
def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments as an object.
    """
    parser = argparse.ArgumentParser()
    # # Check if there is a stdin 
    # if sys.stdin.isatty():
    #     # If not add input dir from arguments
    #     parser.add_argument(
    #         '-id', '--input_dir', type=str, required=True, 
    #         help="Path to the input directory, can also use stdin with no '-id' argument."
    #     )
    parser.add_argument(
        '-m', '--fn_metat', type=str, required=True,
        help="Path to metaT file."
    )
    parser.add_argument(
        '-mt', '--fn_metat_tar', nargs='?', const='',
        help="If the metat file is a tarball, what is the name of the file within the tarball to extract."
    )
    parser.add_argument(
        '-t', '--fn_targets', type=str, required=True,
        help="Path to a file with a list of targets to find."
    )

    parser.add_argument(
        '-k', '--name_key', type=str, nargs='?', const='',
        help="Column name in metat file to be used for dictionary keys."
    )
    parser.add_argument(
        '-vl', '--name_value', type=str, nargs='?', const='',
        help="Column name in metat file to be used for dictionary values."
    )
    parser.add_argument(
        '-c', '--columns_line_number', type=int, nargs='?', const=1,
        help="Line in the file to search for column names. Default is line 1 (counting starts at 1 not 0)."
    )
    parser.add_argument(
        '-ik', '--idx_key', type=int, nargs='?', const=0,
        help="Column name in metat file to be used for dictionary keys."
    )
    parser.add_argument(
        '-ivl', '--idx_value', type=int, nargs='?', const=1,
        help="Column name in metat file to be used for dictionary values."
    )

    parser.add_argument(
        '-o', '--output_fn', type=str, default="", 
        help="Path to the output file. Default is 'output/{fn_metat}-{fn_targets}-dict-{name_key}-{name_value}.json'."
    )
    parser.add_argument(
        '--verbose', action='store_true', 
        help="Enable verbose mode for detailed logging."
    )
    
    args = parser.parse_args()

    # # Check if there is stdin
    # if not sys.stdin.isatty():
    #     # Read all lines from stdin
    #     input_dir = sys.stdin.read() 
    # else: 
    #     # If not add input dir from args  
    #     input_dir = args.input_dir

    return(args)

###################################################################################################
# Main
###################################################################################################

def main():
    """
    Main function to execute script logic.
    """
    # Parse command line arguments
    args = parse_arguments()

    if args.verbose:
        print("Verbose mode enabled.")
        print(f"Input fn: {args.fn_metat}")
        if args.fn_metat_tar:
            print(f"Tarball fn: {args.fn_metat_tar}")
        print(f"Target fn: {args.fn_targets}")
        print(f"Key name: {args.name_key}")
        print(f"Value name: {args.name_value}")
        print(f"Column line number: {args.columns_line_number}")
        if args.idx_key:
            print(f"Key index: {args.idx_key}")
            print(f"Value index: {args.idx_value}")
        print(f"Output fn: {args.output_fn}")
        

    # Read KO list
    with open(args.fn_targets) as f:
        set_search = set(f.read().splitlines())

    # Get key, value indices
    if args.name_key and args.name_value:
        idx_key, idx_value = get_key_value_idx(
            args.fn_metat, 
            args.fn_metat_tar,
            args.name_key.replace('"','').replace("'", ""), 
            args.name_value.replace('"','').replace("'", ""), 
            args.columns_line_number
        )
    else:
        idx_key = args.idx_key
        idx_value = args.idx_value

    # Get dictionary from metat file, mapping key to value
    dict_metat = get_metat_dict(
        args.fn_metat,
        args.fn_metat_tar,
        set_search,
        idx_key, 
        idx_value,   
    )
    
    # Get output filename
    if not args.output_fn:
        out_fn = get_output_fn(args)
    else:
        out_fn = args.output_fn

    # Write dictionary to json
    save_dict(dict_metat, out_fn)



if __name__ == "__main__":
    main()
