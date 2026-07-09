process ARRIBA_FUSION_ALIGNMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fusions_tsv)

    output:
    tuple val(meta), path("*_fusions_alignments.bam"), path("*_fusions_alignments.bam.bai"), emit: fusion_bam
    path "versions.yml"                                                                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_fusion-supporting_alignments.sh ${fusions_tsv} ${bam} extracted

    if [ -e extracted_1.bam ]; then
        samtools merge -O SAM - extracted*.bam
    else
        samtools view -H ${bam}
    fi | \\
    awk '!duplicate[\$0]++' | \\
    samtools view -O "BAM,level=9" --write-index -o ${prefix}_fusions_alignments.bam##idx##${prefix}_fusions_alignments.bam.bai

    # Clean up intermediate files to save disk space
    rm -f extracted*.bam*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}