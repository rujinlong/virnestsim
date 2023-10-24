#!/usr/bin/env nextflow
// Usage: nextflow run rujinlong/virnestsim -profile denglab_hmgu,singularity -resume

nextflow.enable.dsl=2

log.info """\
NF - VirNest simulation pipeline
=========================
result: ${params.outdir}
report: ${params.report}
"""

process CAMISIM {
    input:
    path camisimcfg
    path camisimmeta
    path camisimtsv
    path refseqs

    output:
    path "output.txt"

    script:
    """
    metagenomesimulation.py $camisimcfg
    """
}


workflow {
    // Contigs from Lab isolation and/or NCBI genomes
    camisimcfg_ch = Channel.fromPath(params.camisimcfg)
    camisimmeta_ch = Channel.fromPath(params.camisimmeta)
    camisimtsv_ch = Channel.fromPath(params.camisimtsv)
    refseqs_ch = Channel.fromPath(params.refseqs)
    CAMISIM(camisimcfg_ch, camisimmeta_ch, camisimtsv_ch, refseqs_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
