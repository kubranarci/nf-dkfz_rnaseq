//
// RUN_ARRIBA
//

include { ARRIBA_ARRIBA            } from '../../../modules/nf-core/arriba/arriba/main'
include { ARRIBA_VISUALISATION     } from '../../../modules/nf-core/arriba/visualisation/main'
include { SAMTOOLS_SORT            } from '../../../modules/nf-core/samtools/sort/main'
include { ARRIBA_QUANTIFY_VIRUS    } from '../../../modules/nf-core/arriba/quantify_virus/main'
include { ARRIBA_VIRAL_ALIGNMENTS  } from '../../../modules/nf-core/arriba/viral_alignments/main'
include { ARRIBA_FUSION_ALIGNMENTS } from '../../../modules/nf-core/arriba/fusion_alignments/main'


workflow ARRIBA_WF {
    take:
    reads               // channel: [val(meta), bam]
    fasta               // channel: [val(meta), fasta]
    gencode_gtf_ch      // channel: [val(meta), gtf]
    blacklist           // channel: [file]
    known_fusions       // channel: [file]
    cytobands           // channel: [file]
    protein_domains     // channel: [file]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    SAMTOOLS_SORT(
        reads,
        fasta,
        'bai'
    )

    ARRIBA_ARRIBA (
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai).map{meta, bam, bai -> tuple(meta, bam, bai, [])},
        fasta,
        gencode_gtf_ch,
        blacklist,
        known_fusions,
        cytobands,
        protein_domains
    )
    versions = versions.mix(ARRIBA_ARRIBA.out.versions)
    multiqc_files = multiqc_files.mix(ARRIBA_ARRIBA.out.fusions_fail)

    bams_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai)

    if (params.viralContigs){
        ARRIBA_QUANTIFY_VIRUS(
            bams_ch
            )
        
        ARRIBA_VIRAL_ALIGNMENTS(
            bams_ch
        )

        bams_ch = ARRIBA_VIRAL_ALIGNMENTS.out.viral_bam

    }

    ARRIBA_FUSION_ALIGNMENTS(
        bams_ch.join(ARRIBA_ARRIBA.out.fusions)
    )

    ARRIBA_VISUALISATION(
        bams_ch.join(ARRIBA_ARRIBA.out.fusions),
        gencode_gtf_ch,
        protein_domains.map{file -> tuple([id:"protein_domains"], file)},
        cytobands.map{file -> tuple([id:"cytobands"], file)}
    )
    multiqc_files = multiqc_files.mix(ARRIBA_VISUALISATION.out.pdf)

    emit:
    versions
    multiqc_files

}