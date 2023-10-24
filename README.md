# virnestsim

# Usage

**Download NCBI taxonomy**

```bash
wget https://hmgubox2.helmholtz-munich.de/index.php/s/5BFYRAxnWemw8e6/download/NCBI.tar.gz
```

**Run the pipeline**

```bash
# download the latest version pipeline
nextflow pull rujinlong/virnestsim

# run the pipeline
nextflow run rujinlong/virnestsim \
    -r dev \
    -profile standard,singularity \
    --camisimcfg virnest_contig.ini \
    --camisimmeta metadata_virnest.tsv \
    --camisimtsv genome_to_id_virnest.tsv \
    --refseqs input_virclust \
    --ncbi NCBI.tar.gz
```
