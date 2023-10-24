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

    output:
    path "output.txt"

    script:
    """
    mategenomesimulation.py $camisimcfg
    """
}


workflow {
    // Contigs from Lab isolation and/or NCBI genomes
    camisimcfg_ch = Channel.fromPath(params.camisimcfg)
    CAMISIM(camisimcfg_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
