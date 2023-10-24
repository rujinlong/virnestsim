# virnestsim

# Usage

```bash
# download the latest version pipeline
nextflow pull rujinlong/virnestsim

# run the pipeline
nextflow run rujinlong/virnestsim \
    -r dev \
    -profile standard,singularity \
    --camisimcfg "virnest_contig.ini" \
    --camisimmeta metadata_virnest.tsv \
    --camisimtsv genome_to_id_virnest.tsv \
    --refseqs input_virclust \
    --ncbi NCBI.tar.gz
```
