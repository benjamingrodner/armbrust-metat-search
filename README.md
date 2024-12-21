# Armbrust metaT Search

Tools for searching the Armbrust lab Gradients metaT dataset.

Currently set up for dictionary construction. Search each line of the tables for a set of keys and map to the corresponding value in the line.

## Docker setup

Custom build:

```
docker build --platform linux/amd64 --tag <username>/<image> .
docker push <username>/<image>
```

Otherwise pull from docker://benjamingrodner/get_metat_dicts

## Example usage

Map a set of Kegg orthologies to the corresponding contigs.

```
cd example

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
````