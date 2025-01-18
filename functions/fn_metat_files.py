"""
Package Name: fn.metat_files.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2025-01-17
Description:
    Standard functions to manipulate metat files

Usage:
    import sys
    sys.path.append(<path to repo>/scripts/fn_metat_files.py)
    import fn_metat_files as fmf

"""

import os
import re
import tarfile
import gzip

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
        fspl, ext = os.path.splitext(fn_tar)
        if ext == '.gz':
            ext = os.path.splitext(fspl)[1]
        line = line.decode('utf-8')
    # Separate based on file type
    if ext == '.csv':
        return re.split(r',', line)
    if ext == '.tsv':
        return re.split(r'\t', line)
    if ext == '.tab':
        return re.split(r'\s+', line)
    


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
        # If its a tarball
        tarf = tarfile.open(fn, "r|gz")
        ext = os.path.splitext(fn_tar)[1]
        for t in tarf:
            if fn_tar in t.name:
                f = tarf.extractfile(t)
                # If the files within the tarball are gzipped
                if ext == '.gz':
                    f = gzip.open(f)
                return f
        # If the file failed to read
        tarf = tarfile.open(fn, "r|gz")
        ext = os.path.splitext(fn_tar)[1]
        i = 0
        tnames = []
        for t in tarf:
            tnames.append(t.name)
            i += 1
            if i > 5:
                break
        raise ValueError(
            f"""
            The target subfilename failed to match any files in the tarball. 
            The target subfilename was:
            {fn_tar}
            Here are 5 tar subfilenames from the tarball:
            {tnames}
            The tarball filename was:
            {fn}
            """
        )


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


def get_column_name_idx(fn, fn_tar, name, cln=1):
    """
    Get indices for keys and values in parsed metat file.

    Args:
        fn (str): Path to metat file.
        fn_tar (str): If the metat file is a tarball, what is the name of the file within the tarball to extract.
        name (str): Column name in metat file to be used for dictionary values.

    Returns:
        int: Index of the column name in the parsed line (parsed by ',' or '\s').
    """
    f = open_file(fn, fn_tar)
    for _ in range(cln):
        line = f.readline()
    f.close()
    l = split_line(line, fn, fn_tar)
    for idx, column in enumerate(l):
        # Remove extra quotes
        column = column.replace('"','').replace("'", "")
        if column == name:
            idx_name = idx

    try:
        idx_name
    except NameError:
        raise ValueError(
            f"The provided argument -k {name} did not map to column names in the file.\n\
            Here is the split line (line number {cln} in the file) used to search for column names:\n\
            {l}\n\
            To change the column line, specify -c <columns_line_number> in the arguments (counting starts from 1 not 0)"
        )

    return idx_name