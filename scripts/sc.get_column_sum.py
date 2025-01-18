#!/usr/bin/env python

"""
Script Name: sc.get_column_sum.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2025-01-17
Description:
    Given a metat file and a column
    Sum the values in that column

Usage:
    python sc.get_column_sum.py -m <fn_metat> -vl <name_value> -c <columns_line_number> -ivl <idx_value> -o <output_fn> --verbose

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

# Custom functions
import fn_metat_files as fnf

###################################################################################################
# Functions
###################################################################################################
def save_file(c, fn):
    """
    Write a human readable json file

    Args:
        c (float): Value to write.
        fn (str): Output filename.
    """
    # Make folder if not existing
    od = os.path.split(fn)[0]
    if not os.path.exists(od):
        os.makedirs(od)
        print(f"Made dir: {od}")
    with open(fn, 'w') as f:
        f.write(str(c))
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
    col = args.name_value or args.idx_value
    return(
        'output/' 
        + args.fn_metat + '-' + 
        + '.sum_column_' + col +
        + '.txt'
    )

        
def get_column_sum(
        fn,
        fn_tar,
        idx=0
    ):
    """
    Sum values in the column of a metat file
    SET UP to parse armbrust-metat folder files as of 2025-01-17

    Args:
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
        idx (int): Index of the colum in the parsed line.

    Returns:
        float: Sum of values that can be interpreted as float in the column.
    """
    count = 0
    i = 0
    f = fnf.open_file(fn, fn_tar)
    for line in f:
        l = fnf.split_line(line, fn, fn_tar)
        if len(l) > idx:
            # Get key and remove extra quotes
            value = l[idx].replace('"','').replace("'", "")
            try:
                count += float(value)
                i += 1
            except:
                pass
    f.close()
    if not i:
        raise ValueError(
            f"""
            The target column did not have any values that could be interpreted as float.
            Target column idx: {idx}
            File name: {fn}
            File name tar: {fn_tar}
            """
        )
    return count




    
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
        '-vl', '--name_value', type=str, nargs='?', const='',
        help="Column name in metat file to be used for dictionary values."
    )
    parser.add_argument(
        '-c', '--columns_line_number', type=int, nargs='?', const=1,
        help="Line in the file to search for column names. Default is line 1 (counting starts at 1 not 0)."
    )
    parser.add_argument(
        '-ivl', '--idx_value', type=int, nargs='?', const=1,
        help="Column name in metat file to be used for dictionary values."
    )
    # parser.add_argument(
    #     '-op', '--operation', type=int, nargs='?', const='sum',
    #     help="Operation to perform on the colum. Options are: 'sum','mean','std'"
    # )
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
        print(f"Value name: {args.name_value}")
        print(f"Column line number: {args.columns_line_number}")
        if args.idx_value:
            print(f"Value index: {args.idx_value}")
        # print(f"Operation: {args.operation}")
        print(f"Output fn: {args.output_fn}")
        

    # Get key, value indices
    if args.name_value:
        idx_value = fnf.get_column_name_idx(
            args.fn_metat, 
            args.fn_metat_tar,
            args.name_value.replace('"','').replace("'", ""), 
            args.columns_line_number
        )
    else:
        idx_value = args.idx_value

    # Sum column
    column_sum = get_column_sum(
        fn=args.fn_metat,
        fn_tar=args.fn_metat_tar,
        idx=idx_value,
    )
    
    # Get output filename
    if not args.output_fn:
        out_fn = get_output_fn(args)
    else:
        out_fn = args.output_fn

    # Write dictionary to json
    save_file(column_sum, out_fn)



if __name__ == "__main__":
    main()


