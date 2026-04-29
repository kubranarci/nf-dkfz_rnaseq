process DRAW_FUSIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    tuple val(meta), path(fusions)
    tuple val(meta2),path (annotation)

    output:
    tuple val(meta), path("*.fusions.pdf"), emit: pdf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    draw_fusions.R \\
        --annotation=${annotation} \\
        --fusions=${fusions} \\
        --output=${prefix}.fusions.pdf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' | awk '{print \$2}')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fusions.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' | awk '{print \$2}')
    END_VERSIONS
    """
}