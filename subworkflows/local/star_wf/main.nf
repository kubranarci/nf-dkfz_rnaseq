//
// RUN_STAR
//

include { STAR_ALIGN          } from '../../../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'
include { SAMTOOLS_SORT       } from '../../../modules/nf-core/samtools/sort/main'
include { SAMBAMBA_MARKDUP    } from '../../../modules/nf-core/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT   } from '../../../modules/nf-core/sambamba/flagstat/main'

workflow STAR_WF {
    take:
    samplesheet     // channel: [val(meta), reads ]
    index           // channel: [val(meta), star_index] 
    fasta           // channel: [val(meta), fasta]
    gtf             // channel: [val(meta), gtf]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    // Fixed the space in the || operator
    if (!params.star_index) {

        STAR_GENOMEGENERATE (
            fasta,
            gtf
        )
        index = STAR_GENOMEGENERATE.out.index
    }
    
    STAR_ALIGN (
        samplesheet,
        index,
        gtf,
        true,
        params.seq_platform,
        params.seq_center
    )
    versions = STAR_ALIGN.out.versions
    ch_star_transcripts = STAR_ALIGN.out.bam_transcript
    multiqc_files = multiqc_files.mix(STAR_ALIGN.out.log_final)

    ch_star_unsorted_bam = STAR_ALIGN.out.bam_unsorted
    sorted_bam = STAR_ALIGN.out.bam_sorted_aligned.map { meta, sam ->
        def new_meta = meta + [ chimera: false ]
        return [ new_meta, sam ]
    }
    chimera_sam = STAR_ALIGN.out.sam.map { meta, sam ->
        def new_meta = meta + [ chimera: true ]
        return [ new_meta, sam ]
    }

    SAMTOOLS_SORT(
        chimera_sam,
        fasta,
        "bai"
    )

    SAMBAMBA_MARKDUP(
        SAMTOOLS_SORT.out.bam.mix(sorted_bam)
    )   
    sorted_mkdup_bam = SAMBAMBA_MARKDUP.out.bam.join(SAMBAMBA_MARKDUP.out.bai)
    versions         = versions.mix(SAMBAMBA_MARKDUP.out.versions)

    SAMBAMBA_FLAGSTAT(
        SAMBAMBA_MARKDUP.out.bam
    )
    versions      = versions.mix(SAMBAMBA_FLAGSTAT.out.versions)
    multiqc_files = multiqc_files.mix(SAMBAMBA_FLAGSTAT.out.stats)
    flagstat      = SAMBAMBA_FLAGSTAT.out.stats

    sorted_mkdup_bam.branch { meta, _bam, _bai ->
        chimera: meta.chimera == true
        sorted_bam: meta.chimera == false
    }.set{ch_bam}

    // we only need sorted bam output
    ch_star_sorted_mkdup_bam = ch_bam.sorted_bam

    emit:
    versions
    multiqc_files
    ch_star_sorted_mkdup_bam
    ch_star_unsorted_bam
    ch_star_transcripts
    chimera_sam
    flagstat
    index
}