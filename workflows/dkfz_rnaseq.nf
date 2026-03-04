/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dkfz_rnaseq_pipeline'
include { STAR                   } from '../subworkflows/local/star/main'
include { FINGERPRINTING         } from '../modules/local/fingerprinting/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DKFZ_RNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    // create reference channels
    fasta_ch = params.fasta ? channel.fromPath(params.fasta, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    fai_ch   = params.fai   ? channel.fromPath(params.fai, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    star_ch  = params.star  ? channel.fromPath(params.star, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    gtf_ch   = params.gtf   ? channel.fromPath(params.gtf, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    fingerprinting_sites_ch = params.fingerprinting_sites ? channel.fromPath(params.fingerprinting_sites, checkIfExists: true).collect() : channel.empty()

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // add feature: skip_alignmnet
    // 
    if (params.skip_tools.contains("star")) {
        log.warn "Skipping STAR alignment as requested with --skip_tools. Downstream steps that depend on STAR output will likely fail, so use with caution."
    } else {

        STAR (
            ch_samplesheet,
            star_ch,
            fasta_ch,
            gtf_ch
        )
        ch_multiqc_files = ch_multiqc_files.mix(STAR.out.multiqc_files.map{ meta, files -> files })
        ch_versions = ch_versions.mix(STAR.out.versions)

        if (params.skip_tools.contains("fingerprinting")) {
            log.warn "Skipping fingerprinting as requested with --skip_tools. Downstream steps that depend on Kallisto output will likely fail, so use with caution."
        }else {
            FINGERPRINTING(
                STAR.out.ch_star_sorted_mkdup_bam,
                fingerprinting_sites_ch
            )
            ch_versions = ch_versions.mix(FINGERPRINTING.out.versions)
        }

        if (params.skip_tools.contains("rnaseqc")){
            log.warn "Skipping RNASeQC as requested with --skip_tools. Downstream steps that depend on RNASeQC output will likely fail, so use with caution."
        }else{
            //RNASEQC()
        }
    }




    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'dkfz_rnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
