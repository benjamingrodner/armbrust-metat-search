#!/bin/bash

# This runs the Snakefile pipeline for a table of files
# The config.example.yaml file points to the table (file_table.example.csv)
# using a docker image 

apptainer run \
    --no-home \
    docker://benjamingrodner/get_metat_dicts \
    Snakefile \
        --configfile config.example.yaml \
        -j 1 \
        -p