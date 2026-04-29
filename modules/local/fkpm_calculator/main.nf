process FKPM_CALCULATOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::perl=5.32.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/perl:5.26.2'
        : 'quay.io/biocontainers/perl:5.26.2'}"

    input:
    tuple val(meta), path(featurecounts)
    tuple val(meta2), path(summaries)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(gtf_exclude)

    output:
    tuple val(meta), path("*.tsv")         ,emit: fpkm_tpm
    tuple val(meta), path("*_raw.tgz")     , emit: raw_counts_archive
    path  "versions.yml"                   , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # 1. Safely find the script and export the Perl library path
    SCRIPT_PATH=\$(command -v featureCounts_2_FpkmTpm.pl || true)
    
    if [ -z "\$SCRIPT_PATH" ]; then
        echo "ERROR: featureCounts_2_FpkmTpm.pl not found in PATH."
        echo "Did you put it in the bin/ directory and run 'chmod +x'?"
        exit 1
    fi
    
    # Native bash way to get the directory path without using 'dirname'
    export PERL5LIB="\${SCRIPT_PATH%/*}":\${PERL5LIB:-}

    # 2. Run the custom DKFZ Perl script
    featureCounts_2_FpkmTpm.pl \\
        $featurecounts \\
        $gtf \\
        $gtf_exclude \\
        > ${prefix}.tsv

    # 3. Cleanup: Move raw counts to a folder and compress them
    mkdir ${prefix}_raw
    mv $featurecounts ${prefix}_raw/
    mv $summaries ${prefix}_raw/

    # Create the archive without the z option and compress it manually
    tar -cvf ${prefix}_raw.tar ${prefix}_raw/
    gzip ${prefix}_raw.tar
    mv ${prefix}_raw.tar.gz ${prefix}_raw.tgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.tsv
    mkdir ${prefix}_raw
    tar -cvf ${prefix}_raw.tar ${prefix}_raw/
    gzip ${prefix}_raw.tar
    mv ${prefix}_raw.tar.gz ${prefix}_raw.tgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -v | grep 'version' | sed 's/.*(v//; s/).*//' )
    END_VERSIONS
    """  
}