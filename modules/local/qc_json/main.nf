process QC_JSON {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::perl=5.32.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/perl:5.26.2'
        : 'quay.io/biocontainers/perl:5.26.2'}"

    input:
    tuple val(meta), path(rnaseqc_metrics),path(flagstat)

    output:
    tuple val(meta), path("*_qualitycontrol.json"), emit: json
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    qcJson.pl \\
        ${flagstat} \\
        ${rnaseqc_metrics} \\
        > ${prefix}_qualitycontrol.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_qualitycontrol.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """
}