
mkdir -p jupyter_home

apptainer run \
    --home jupyter_home:/app\
    --bind /mnt/nfs/projects/armbrust-metat \
    --bind /scratch/bgrodner \
    docker://benjamingrodner/metat-jupyter \
        jupyter notebook \
            --ip=0.0.0.0 \
            --port=8889 \
            --no-browser \
            --NotebookApp.allow_origin='*'