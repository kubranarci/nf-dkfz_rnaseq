process ARRIBA_VIRAL_ALIGNMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_virus_alignments.bam"), path("*_virus_alignments.bam.bai"), emit: viral_bam
    path "versions.yml"                                                                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export VIRAL_CONTIGS='${params.viralContigs}'

    # Pre check to see if any viral reads exist
    HAS_VIRUS=\$(samtools idxstats ${bam} | awk -v vc="\$VIRAL_CONTIGS" '\$3>0 && match(\$1,vc) {print \$1}' | head -n 1)

    if [ -z "\$HAS_VIRUS" ]; then
        # No viruses found. Generate an empty BAM containing only the header.
        samtools view -H ${bam} | samtools view -O "BAM,level=9" --write-index -o ${prefix}_virus_alignments.bam##idx##${prefix}_virus_alignments.bam.bai
    else
        # Viruses found. Extract them safely.
        # Using awk OFS instead of tabs in the print statement prevents Groovy string interpolation errors
        samtools idxstats ${bam} | \\
        awk -v OFS='\\t' -v contigs="\$VIRAL_CONTIGS" '\$3 > 0 && match(\$1, contigs) {print \$1, 0, \$2}' | \\
        samtools view -F 2048 -M -L /dev/stdin -O "BAM,level=9" --write-index -o ${prefix}_virus_alignments.bam##idx##${prefix}_virus_alignments.bam.bai ${bam}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}