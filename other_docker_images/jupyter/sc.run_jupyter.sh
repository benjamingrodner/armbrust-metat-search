apptainer run \
    --no-home \
    --bind /mnt/nfs/projects/armbrust-metat \
    --bind /scratch/bgrodner \
    docker://benjamingrodner/metat-jupyter \
        jupyter notebook \
            --ip=0.0.0.0 \
            --port=8889 \
            --no-browser \
            --NotebookApp.allow_origin='*'