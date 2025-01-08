#!/usr/bin/env python

"""
Script Name: sc_dict_values_to_list.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2024-01-03
Description:
    Given a file with a dictionary of lists
    Write a file with the merged lists

Usage:
    python sc_dict_values_to_list.py 

"""

###################################################################################################
# Imports
###################################################################################################

import argparse
import json
import os

###################################################################################################
# Functions
###################################################################################################


def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments as an object.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', '--fn_dict_lists', type=str, required=True,
        help="Path to dict of lists file."
    )

    parser.add_argument(
        '-o', '--fn_out', type=str, default="", 
        help="Path to the output file. Default is 'output/{fn_dict_lists}-list.txt'."
    )
    parser.add_argument(
        '--verbose', action='store_true', 
        help="Enable verbose mode for detailed logging."
    )
    
    args = parser.parse_args()

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

    # Set up output file
    if not args.fn_out:
        fn_out = f'output/{args.fn_dict_lists}-list.txt'
    else:
        fn_out = args.fn_out
    dir_out = os.path.split(fn_out)[0]
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    if args.verbose:
        print("Verbose mode enabled.")
        print(f"Input fn: {args.fn_dict_lists}")
        print(f"Output fn: {fn_out}")


    # Load dict
    with open(args.fn_dict_lists, 'r') as f:
        dict_lists = json.load(f)
    
    # Merge and save
    with open(fn_out, 'w') as fo:
        for _, l in dict_lists.items():
            for i in l:
                fo.write(f'{i}\n')


if __name__ == "__main__":

    main()
