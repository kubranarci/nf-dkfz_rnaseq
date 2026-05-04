//
// RUN_FEATURECOUNTS
//

include { SUBREAD_FEATURECOUNTS  } from '../../../modules/nf-core/subread/featurecounts/main'
include { FKPM_CALCULATOR        } from '../../../modules/local/fkpm_calculator/main'
include { FKPM_CALCULATOR_DEXSEQ } from '../../../modules/local/fkpm_calculator_dexseq/main'
 
workflow FEATURECOUNTS_WF {
    take:
    mkdup_map           // channel: [val(meta), bam, bai ]
    gencode_gtf_ch      // channel: [val(meta), gtf]
    gencode_dexseq_ch   // channel: [val(meta), gtf]
    gencode_exclude_ch  // channel: [val(meta), gtf]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    ch_featurecounts_scatter = mkdup_map.flatMap { meta, bam, bai ->
                return [
                    [ meta + [id: "${meta.id}_s0", strandedness: 0], bam, bai ],
                    [ meta + [id: "${meta.id}_s1", strandedness: 1], bam, bai ],
                    [ meta + [id: "${meta.id}_s2", strandedness: 2], bam, bai ]
                ]
            }  

    ch_gtf_to_use = gencode_dexseq_ch
            .concat(gencode_gtf_ch)
            .first()
    

    SUBREAD_FEATURECOUNTS(
        ch_featurecounts_scatter,
        ch_gtf_to_use            
    )
    versions = versions.mix(SUBREAD_FEATURECOUNTS.out.versions)
    multiqc_files = multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary)
    

    ch_featurecounts_gathered = SUBREAD_FEATURECOUNTS.out.counts
                .map { meta, counts ->
                    def original_id = meta.id.replaceAll(/_s[012]$/, "")
                    def new_meta = meta + [id: original_id]
                    new_meta.remove('strandedness')
                    
                    return [ new_meta, counts ]
                }
                .groupTuple(by: 0) 

    ch_summaries_gathered = SUBREAD_FEATURECOUNTS.out.summary
                .map { meta, summary ->
                    def original_id = meta.id.replaceAll(/_s[012]$/, "")
                    def new_meta = meta + [id: original_id]
                    new_meta.remove('strandedness')
                    
                    return [ new_meta, summary ]
                }
                .groupTuple(by: 0) 


    ch_run_dexseq_flag = gencode_dexseq_ch
        .map { true }
        .ifEmpty( false )

    ch_routed_counts = ch_featurecounts_gathered
        .combine(ch_run_dexseq_flag)
        .branch { _meta, _counts, is_dexseq ->
            dexseq: is_dexseq == true
            standard: is_dexseq == false
        }

    ch_routed_summaries = ch_summaries_gathered
        .combine(ch_run_dexseq_flag)
        .branch { _meta, _summary, is_dexseq ->
            dexseq: is_dexseq == true
            standard: is_dexseq == false
        }

    
    FKPM_CALCULATOR_DEXSEQ(
        ch_routed_counts.dexseq.map { meta, counts, _flag -> [meta, counts] },
        ch_routed_summaries.dexseq.map { meta, summary, _flag -> [meta, summary] },
        gencode_gtf_ch, 
        gencode_exclude_ch
    )
    versions = versions.mix(FKPM_CALCULATOR_DEXSEQ.out.versions)
    multiqc_files = multiqc_files.mix(FKPM_CALCULATOR_DEXSEQ.out.fpkm_tpm)



    FKPM_CALCULATOR(
        ch_routed_counts.standard.map { meta, counts, _flag -> [meta, counts] },
        ch_routed_summaries.standard.map { meta, summary, _flag -> [meta, summary] },
        gencode_gtf_ch,
        gencode_exclude_ch
    )
    versions = versions.mix(FKPM_CALCULATOR.out.versions)
    multiqc_files = multiqc_files.mix(FKPM_CALCULATOR.out.fpkm_tpm)

    emit:
    versions
    multiqc_files

}