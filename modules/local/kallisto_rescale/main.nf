process KALLISTO_RESCALE {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::perl=5.32.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/perl:5.26.2'
        : 'quay.io/biocontainers/perl:5.26.2'}"

    input:
    tuple val(meta), path(abundance)
    tuple val(meta2), path(gtf_exclude)

    output:
    tuple val(meta), path("*.tsv")         , emit: rescaled
    path  "versions.yml"                   , emit: versions

    script:
    """
    kallistoRescaleTPM.pl \\
        $abundance/abundance.tsv \\
        $gtf_exclude \\
        > abundance.rescaled.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """

    stub:
    """
    touch abundance.rescaled.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """  
}