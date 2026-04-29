//
// RUN_QC
//

include { FINGERPRINTING    } from '../../../modules/local/fingerprinting/main'
include { RNASEQC           } from '../../../modules/local/rnaseqc/main'
include { QUALIMAP_RNASEQ   } from '../../../modules/nf-core/qualimap/rnaseq/main'
include { QC_JSON           } from '../../../modules/local/qc_json/main'

workflow QC_WF {
    take:
    ch_star_sorted_mkdup_bam     // channel: [val(meta), bam, bai ]
    ch_star_unsorted_bam         // channel: [val(meta), bam, bai] 
    fasta                        // channel: [val(meta), fasta]
    fai                          // channel: [val(meta), fai]
    fingerprinting               // channel: bed
    gtf                          // channel: [val(meta), gtf]
    flagstat                     // channel: [val(meta), txt]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    if (params.skip_tools.contains("fingerprinting")) {
        log.warn "Skipping fingerprinting as requested with --skip_tools. Downstream steps that depend on fingerprinting output will likely fail, so use with caution."
    }else {
        FINGERPRINTING(
            ch_star_sorted_mkdup_bam,
            fingerprinting
        )
        versions = versions.mix(FINGERPRINTING.out.versions)
    }

    if (params.skip_tools.contains("rnaseqc")){
        log.warn "Skipping RNASeQC as requested with --skip_tools. Downstream steps that depend on RNASeQC output will likely fail, so use with caution."
    }else{
        RNASEQC(
            ch_star_sorted_mkdup_bam,
            gtf,
            fasta,
            fai
        )
    }

    if (params.skip_tools.contains("qualimap")){
        log.warn "Skipping QUALIMAP as requested with --skip_tools. Downstream steps that depend on QUALIMAP output will likely fail, so use with caution."
    }else{
        QUALIMAP_RNASEQ(
            ch_star_unsorted_bam,
            gtf
        )
        versions = versions.mix(QUALIMAP_RNASEQ.out.versions)  
    }

    if (params.skip_tools.contains("qcjson")){
        log.warn "Skipping QCJSON as requested with --skip_tools. Downstream steps that depend on QCJSON output will likely fail, so use with caution."
    }else{
        QC_JSON(
            RNASEQC.out.metrics.join(flagstat)
        )
        versions = versions.mix(QC_JSON.out.versions)  
    }


    emit:
    versions
    multiqc_files

}