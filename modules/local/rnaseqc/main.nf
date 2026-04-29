process RNASEQC {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::rna-seqc=2.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/42/429086de4250328c3aa597731c9088c27dd5e536f6a199ed33218f6d62c408f0/data' :
        'community.wave.seqera.io/library/rna-seqc:2.4.2--1db386cefd58a1f8' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.gene_reads.gct")   , emit: counts
    tuple val(meta), path("*.gene_tpm.gct")     , emit: tpm
    tuple val(meta), path("*.metrics.tsv")      , emit: metrics
    tuple val(meta), path("*.exon_reads.gct")   , optional: true, emit: exon_counts
    tuple val(meta), path("*.coverage.tsv")     , optional: true, emit: coverage
    tuple val("${task.process}"), val('rnaseqc'), eval("rnaseqc --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_rnaseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-u" : ""
    """
    rnaseqc \\
        $gtf \\
        $bam \\
        . \\
        $args \\
        $paired \\
        --fasta $fasta \\
        --sample $prefix
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_reads.gct
    touch ${prefix}.gene_tpm.gct
    touch ${prefix}.metrics.tsv
    touch ${prefix}.exon_reads.gct
    touch ${prefix}.coverage.tsv
    """             
}