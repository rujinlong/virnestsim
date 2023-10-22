#!/usr/bin/env nextflow
// Usage: nextflow run rujinlong/virnestsim -profile denglab_hmgu,singularity -resume

nextflow.enable.dsl=2

log.info """\
NF - VirNest simulation pipeline
=========================
result: ${params.outdir}
report: ${params.report}
"""

process copyFile {

    input:
    path inputFile

    output:
    path "output.txt"

    script:
    """
    cat $inputFile > output.txt
    """
}


workflow plass {
    // Viral-tagging assembly (plass protein assembly + plass guided_nuclassemble)
    take:
    raw_reads_ch

    main:
    FASTQC(raw_reads_ch)
    CLEAN_READS(raw_reads_ch)

    ASSEMBLY_PLASS(CLEAN_READS.out.clean_reads_paired_ch)
    ASSEMBLY_PENGUIN(CLEAN_READS.out.clean_reads_paired_ch)

    if (params.mode == "fastqc") {
        MULTIQC(FASTQC.out.fastqc_results_ch.collect())
    }
    else {
        MULTIQC(FASTQC.out.fastqc_results_ch.flatten().concat(CLEAN_READS.out.fastp_json_ch.flatten()).collect())
    }

    emit:
    FAA_PLASS = ASSEMBLY_PLASS.out.faa_plass_ch
    CONTIG_PLASS = ASSEMBLY_PENGUIN.out.penguin_ch
}

workflow {
    // Contigs from Lab isolation and/or NCBI genomes
    contigs_ch = Channel.fromPath(params.contigs)
    // Proteins from NCBI genomes (can be empty?)
    prot_custom = Channel.fromPath(params.prot_custom)
    gene_custom = Channel.fromPath(params.gene_custom)

    // PLASS prot and hybrid assembly
    if (params.reads != 'null') {
        raw_reads_ch = Channel.fromFilePairs(params.reads)
        plass(raw_reads_ch)
        // DVF(plass.out.CONTIG_PLASS)
        GENOMAD(plass.out.CONTIG_PLASS)
        // contigs_ch2 = contigs_ch.mix(DVF.out.viral_contigs_ch).collect()
        contigs_ch3 = contigs_ch.mix(GENOMAD.out.viral_contigs_ch).collect()
    } else {
        // contigs_ch2 = contigs_ch
        contigs_ch3 = contigs_ch
    }

    // isolated genomes: Lab genomes (and/or) NCBI genomes
    // PROT_PRED(contigs_ch2)
    PROT_PRED(contigs_ch3)
    
    // Merge proteins
    if (params.reads == 'null') {
        prot_ch = PROT_PRED.out.prot_pred_faa.mix(prot_custom).collect()
    } else {
        prot_ch = PROT_PRED.out.prot_pred_faa.mix(plass.out.FAA_PLASS).mix(prot_custom).collect()
    }
    // Merge all proteins -> PCs
    PROT_CLUSTER(prot_ch)
    PROT_ANNO_PHROG(PROT_CLUSTER.out.prot_cluster_repseqs)

    // Merge genes
    gene_ch = PROT_PRED.out.prot_pred_fna.mix(gene_custom).collect()
    // Merge all genes -> GCs
    GENE_CLUSTER(gene_ch)

    // Viral clusters
    VIR_CLUST_ANI(PROT_PRED.out.contigs_isolates_ch)
    VCONTACT2(VIR_CLUST_ANI.out.vcrep_long_ch)
    MERGE_GENOME_CLUSTERS(VIR_CLUST_ANI.out.species_cluster_ch, VCONTACT2.out.vcontact2_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
