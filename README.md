# Armbrust metaT Search

Tools for searching the Armbrust lab Gradients metaT dataset.

Currently set up for dictionary construction. Search each line of the metaT tables for a set of keys and map hits to the corresponding value in the line. Result is a .json file with dictionary mapping each key -> list of values. 

"Snakefile" - given a list of Kegg orthologies, find all contigs annotated with those KOs.

"Snakefile_find_taxon" - given a list of taxa, find all contigs annotated as those taxa.

"Snakefile_all_taxa_estcounts" - For all taxa, sum all estcounts for all contigs in each taxon.

## Docker setup

Custom build:

```
cd other_docker_images/builder
docker build --platform <> --tag <username>/<image-builder> .
cd ../..
```
Now edit the file "Docker" so that the builder image loads from your build.
```
docker build --platform <> --tag <username>/<image> .
docker push <username>/<image>
```

Otherwise pull from docker://benjamingrodner/get_metat_dicts, which is built with --platform linux/amd64

## Example usage

Map a set of Kegg orthologies to the corresponding contigs.

### Search a single file

```
$ cd example

$ apptainer run \
    --no-home \
    docker://benjamingrodner/get_metat_dicts \
    sc_get_metat_dict.py \
        -m test_table.G1NS.S11C1_3um.kofam2021.subset.csv \
        -t iron_KOs.example.txt \
        -k query_name \
        -vl target_name  \
        -o dicts_iron_KO_contig_example \
        --verbose
```

```
options:
  -h, --help            show this help message and exit

  -m , --fn_metat 
                        Path to the metaT file.
  -t , --fn_targets 
                        Path to a file with a list of targets to find.
  -k , --name_key 
                        Column name in metat file to be used for dictionary keys.
  -vl , --name_value 
                        Column name in metat file to be used for dictionary values.
  -c , --columns_line_number 
                        Line in the file to search for column names. Default is line 1 (counting starts
                        at 1 not 0).
  -o , --output_fn 
                        Path to the output file. Default is
                        'output/{fn_metat}-{fn_targets}-dict-{name_key}-{name_value}.json'.

  --verbose             Enable verbose mode for detailed logging.
```

### Use snakemake to search many files

Requires a configuration file as .yaml or .json, and an input table as .csv specifying the files to search and a very specific set of parameters required to run the pipeline. Examples of these files are shown in the 'example' folder. Example data is also included to test the pipeline.


```
$ cd example

$ apptainer run \
    --no-home \
    -S $HOME
    docker://benjamingrodner/get_metat_dicts \
    Snakefile \
        --configfile config.example.yaml \
        -j 1 \
        -p
```
