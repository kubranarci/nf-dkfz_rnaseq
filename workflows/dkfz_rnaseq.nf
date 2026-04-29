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
include { STAR_WF                } from '../subworkflows/local/star_wf/main'
include { KALLISTO_WF            } from '../subworkflows/local/kallisto_wf/main'
include { ARRIBA_WF              } from '../subworkflows/local/arriba_wf/main'
include { QC_WF                  } from '../subworkflows/local/qc_wf/main'
include { FEATURECOUNTS_WF as FEATURECOUNTS_PLAIN   } from '../subworkflows/local/featurecounts_wf/main'
include { FEATURECOUNTS_WF as  FEATURECOUNTS_DEXSEQ } from '../subworkflows/local/featurecounts_wf/main'

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
    gencode_gtf_ch          = params.gencode_gtf   ? channel.fromPath(params.gencode_gtf, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    gencode_exclude_ch      = params.gencode_exclude   ? channel.fromPath(params.gencode_exclude, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    gencode_dexseq_ch       = params.gencode_dexseq   ? channel.fromPath(params.gencode_dexseq, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()
    fingerprinting_sites_ch = params.fingerprinting_sites ? channel.fromPath(params.fingerprinting_sites, checkIfExists: true).collect() : channel.empty()
    kallisto_ch             = params.kallisto_index  ? channel.fromPath(params.kallisto_index, checkIfExists: true).map { file -> tuple([id: file.getSimpleName()], file) }.collect() : channel.empty()

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // add feature: skip_alignmnet
    // TODO: start from bam files skipping star
    // 
    if (params.skip_tools.contains("star")) {
        log.warn "Skipping STAR alignment as requested with --skip_tools. Downstream steps that depend on STAR output will likely fail, so use with caution."
    } else {
        STAR_WF (
            ch_samplesheet,
            star_ch,
            fasta_ch,
            gtf_ch
        )
        ch_multiqc_files = ch_multiqc_files.mix(STAR_WF.out.multiqc_files.map{ _meta, files -> files })
        ch_versions = ch_versions.mix(STAR_WF.out.versions)
    }

    if (params.skip_tools.contains("qc")) {
        log.warn "Skipping STAR alignment as requested with --skip_tools. Downstream steps that depend on STAR output will likely fail, so use with caution."
    } else {
        QC_WF(
            STAR_WF.out.ch_star_sorted_mkdup_bam,
            STAR_WF.out.ch_star_unsorted_bam,
            fasta_ch,
            fai_ch,
            fingerprinting_sites_ch,
            gencode_gtf_ch,
            STAR_WF.out.flagstat
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_WF.out.multiqc_files.map{ _meta, files -> files })
        ch_versions = ch_versions.mix(QC_WF.out.versions)
    }

    if (params.skip_tools.contains("featurecounts")){
        log.warn "Skipping featurecount as requested with --skip_tools. Downstream steps that depend on featurecount output will likely fail, so use with caution."
    }else{
        FEATURECOUNTS_PLAIN(
            STAR_WF.out.ch_star_sorted_mkdup_bam,
            gencode_gtf_ch,
            channel.empty(),    
            gencode_exclude_ch
        )
        ch_versions = ch_versions.mix(FEATURECOUNTS_PLAIN.out.versions)
           
    }

    if (params.skip_tools.contains("dexseq")){
        log.warn "Skipping dexseq as requested with --skip_tools. Downstream steps that depend on featurecounts_dexseq output will likely fail, so use with caution."
    }else{
        FEATURECOUNTS_DEXSEQ(
            STAR_WF.out.ch_star_sorted_mkdup_bam,
            gencode_gtf_ch,
            gencode_dexseq_ch,
            gencode_exclude_ch
        )
        ch_versions = ch_versions.mix(FEATURECOUNTS_DEXSEQ.out.versions)

    }

    if (params.skip_tools.contains("kallisto")) {
        log.warn "Skipping KALLISTO alignment as requested with --skip_tools. Downstream steps that depend on KALLISTO output will likely fail, so use with caution."
    } 
    else 
    {
        // Branch the trimmed fastq channel to isolate paired-end samples
        // running kallisto with single end reads is still possible if fragment_length and fragment_length_sd given!!
        ch_reads_for_kallisto = ch_samplesheet.branch { meta, _reads ->
            single_end: meta.single_end
            paired_end: !meta.single_end
        }
        
        KALLISTO_WF(
            ch_reads_for_kallisto.paired_end,
            kallisto_ch,
            fasta_ch,
            gencode_exclude_ch
        )
        ch_versions = ch_versions.mix(KALLISTO_WF.out.versions)

    }

    if (params.skip_tools.contains("arriba")){
        log.warn "Skipping Arriba as requested with --skip_tools. Downstream steps that depend on Arriba output will likely fail, so use with caution."
    }
    else
    {
        STAR_WF.out.chimera_sam.map { meta, sam ->
            def new_meta = meta.clone()
            new_meta.remove('chimera')
            return [ new_meta, sam ]
        }.set{ch_star_chimera}

        STAR_WF.out.ch_star_sorted_mkdup_bam.map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.remove('chimera')
            return [ new_meta, bam, bai ]
        }.set{ch_star_sorted_mkdup_bam}

        ARRIBA_WF(
            ch_star_sorted_mkdup_bam.join(ch_star_chimera),
            fasta_ch,
            gencode_gtf_ch,
            [],[]
        )
        ch_versions = ch_versions.mix(ARRIBA_WF.out.versions)

    }


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
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
