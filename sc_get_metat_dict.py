#!/usr/bin/env python

"""
Script Name: sc_get_ko_contig_dicts.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2024-12-19
Description:
    Given a file with a list of Kegg orthologies
    Extract contigs from the gradients dataset with those kegg orthologies

Usage:
    python sc_get_ko_contig_dicts.py -m <fn_metat> -t <fn_targets> -k <name_key> -v <name_value> -c <columns_line_number> -o <output_dir> --verbose

"""

###################################################################################################
# Imports
###################################################################################################

from collections import defaultdict
from time import time
import argparse
import gzip
import json
import sys
import os
import re

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
    bns = []
    for fn_full in [args.fn_metat, args.fn_targets]:
        fn = os.path.split(fn_full)[1]
        bns.append(os.path.splitext(fn)[0])
    return(
        args.output_dir + '/' 
        + bns[0] + '.' + bns[1] 
        + '.dict_' + args.name_key + '_' + args.name_value 
        + '.json'
    )


def split_line(line):
    """
    Split line string read from file

    Args:
        line (str): Unparsed line from the file.
    Returns:
        list: Split line.
    """
    return re.split(r'\s|,', line)


def add_to_dict(line, set_search, idx_key, idx_value, dict_metat):
    """
    Extract key and value from a line and add to the dictionary

    Args:
        line (str): Unparsed line from the file.
        set_search (set): Set of keys to add to the dictionary.
        idx_key (int): Index of the key in the parsed line (parsed by ',' or '\s').
        idx_value (int): Index of the mapped value in the parsed line (parsed by ',' or '\s').
        dict_metat (defaultdict(list)): Existing dictionary.

    Returns:
        defaultdict(list): Updated dictionary.
    """
    # Separate the line by spaces (\s) or commas
    l = split_line(line)
    # Get key and remove extra quotes
    key = l[idx_key].replace('"','').replace("'", "")
    # Search target keys
    if key in set_search:
        # Add value to dict
        value = l[idx_value].replace('"','').replace("'", "")
        dict_metat[key].append(value)
    return dict_metat


def open_file(fn):
    """
    Open text file if gzipped, or not

    Args:
        fn (str): Path to metat file.

    Returns:
        _io.TextIOWrapper: Metat file.
    """
    # Get file type
    ext = os.path.splitext(fn)[1]
    # Unzip to read
    if ext == '.gz':
        return gzip.open(fn, 'rt')
    else:
        return open(fn, 'r')


def get_metat_dict(fn, set_search, idx_key, idx_value):
    """
    Collect all contigs with given Kegg Orthologies
    SET UP to parse armbrust-metat folder files as of 2024-12-19

    Args:
        fn (str): Path to metat file.
        set_search (set): Set of keys to add to the dictionary.
        idx_key (int): Index of the key in the parsed line (parsed by ',' or '\s').
        idx_value (int): Index of the mapped value in the parsed line (parsed by ',' or '\s').

    Returns:
        defaultdict(list): Map of metat values key -> values.
    """
    dict_metat = defaultdict(list)
    f = open_file(fn)
    for line in f:
        dict_metat = add_to_dict(
            line, 
            set_search, idx_key, idx_value,
            dict_metat
        )
    f.close()
    return dict_metat


def get_key_value_idx(fn, name_key, name_value, cln=1):
    """
    Get indices for keys and values in parsed metat file.

    Args:
        fn (str): Path to metat file.
        names_key (str): Column name in metat file to be used for dictionary keys.
        names_value (str): Column name in metat file to be used for dictionary values.

    Returns:
        int: Index of the key in the parsed line (parsed by ',' or '\s').
        int: Index of the mapped value in the parsed line (parsed by ',' or '\s').
    """
    f = open_file(fn)
    for _ in range(cln):
        line = f.readline()
    l = split_line(line)
    for idx, column in enumerate(l):
        # Remove extra quotes
        column = column.replace('"','').replace("'", "")
        if column == name_key:
            idx_key = idx
        elif column == name_value:
            idx_value = idx
    try:
        idx_key, idx_value
    except NameError:
        raise ValueError(
            "The provided arguments -k {name_key} and -vl {name_value} did not map to column names in the file.\n\
            Here is the split line (line number {cln} in the file) used to search for column names:\n\
            {l}\n\
            To change the column line specify -c <columns_line_number> in the arguments (counting starts from 1 not 0)"
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
        help="Path to a file with a list of KOs to find."
    )
    parser.add_argument(
        '-t', '--fn_targets', type=str, required=True,
        help="Path to a file with a list of KOs to find."
    )
    parser.add_argument(
        '-k', '--name_key', type=str, required=True,
        help="Column name in metat file to be used for dictionary keys."
    )
    parser.add_argument(
        '-vl', '--name_value', type=str, required=True,
        help="Column name in metat file to be used for dictionary values."
    )
    parser.add_argument(
        '-c', '--columns_line_number', type=int, default=1,
        help="Line in the file to search for column names. Default is line 1 (counting starts at 1 not 0)."
    )
    parser.add_argument(
        '-o', '--output_dir', type=str, default="metat_dicts", 
        help="Path to the output directory. Default is 'KO_contig_dicts'."
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
        print(f"Target fn: {args.fn_targets}")
        print(f"Key name: {args.name_key}")
        print(f"Value name: {args.name_value}")
        print(f"Column line number: {args.columns_line_number}")
        print(f"Output dir: {args.output_dir}")

    # Read KO list
    with open(args.fn_targets) as f:
        set_search = set(f.read().splitlines())
    # Get key, value indices
    idx_key, idx_value = get_key_value_idx(
        args.fn_metat, 
        args.name_key, 
        args.name_value, 
        args.columns_line_number
    )

    # Get dictionary from metat file, mapping key to value
    dict_metat = get_metat_dict(
        args.fn_metat,
        set_search,
        idx_key, 
        idx_value,   
    )
    
    # Get output filename
    out_fn = get_output_fn(args)
    # Write dictionary to json
    save_dict(dict_metat, out_fn)



if __name__ == "__main__":
    main()
