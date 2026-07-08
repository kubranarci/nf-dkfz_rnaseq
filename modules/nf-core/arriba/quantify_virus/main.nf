process ARRIBA_QUANTIFY_VIRUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    tuple val(meta),  path(bam), path(bai)

    output:
    tuple val(meta), path("*.tsv")          , emit: virus_expression
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export VIRAL_CONTIGS='${params.viralContigs}'

    HAS_VIRUS=\$(samtools idxstats ${bam} | awk -v vc="\$VIRAL_CONTIGS" '\$3>0 && match(\$1,vc) {print \$1}' | head -n 1)

    # 3. Route the logic safely
    if [ -z "\$HAS_VIRUS" ]; then
        # Biological normal: No viral reads found. 
        # Create the empty output file with headers so Nextflow succeeds.
        echo -e "VIRUS\\tGENOME_SIZE\\tCOVERED_BASES\\tCOVERED_GENOME_FRACTION\\tHIGH_QUALITY_ALIGNMENTS\\tRPKM" > ${prefix}.virus_expression.tsv
    else
        # Virus found! Safe to run the script because the region list won't be empty.
        quantify_virus_expression.sh \\
            ${bam} \\
            ${prefix}.virus_expression.tsv
    fi

    # 4. Standard versions output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.virus_expression.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
