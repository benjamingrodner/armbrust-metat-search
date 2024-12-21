#!/bin/bash

# This runs the sc_get_metat_dict.py script for a test sample
# using a docker image

apptainer run \
    --no-home \
    docker://benjamingrodner/get_metat_dicts \
    sc_get_metat_dict.py \
        -m test_table.G1NS.S11C1_3um.kofam2021.subset.csv \
        -t iron_KOs.example.txt \
        -k query_name \
        -vl target_name  \
        -o dicts_iron_KO_contig_example \
        --verbose