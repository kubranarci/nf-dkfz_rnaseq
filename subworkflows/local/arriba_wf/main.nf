//
// RUN_ARRIBA
//

include { ARRIBA_ARRIBA  } from '../../../modules/nf-core/arriba/arriba/main'
include { DRAW_FUSIONS   } from '../../../modules/local/draw_fusions/main'


workflow ARRIBA_WF {
    take:
    reads               // channel [val(meta), star_mkdup_bam, bai, chimera_sam]
    fasta               // channel: [val(meta), fasta]
    gencode_gtf_ch      // channel: [val(meta), gtf]
    known_fusions       // channel: [txt]
    blacklist           // channel: [txt]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    // check and learn optional files
    ARRIBA_ARRIBA (
        reads,
        fasta,
        gencode_gtf_ch,
        [],[],[],[]
    )
    versions = versions.mix(ARRIBA_ARRIBA.out.versions)

    DRAW_FUSIONS(
        ARRIBA_ARRIBA.out.fusions,
        gencode_gtf_ch
    )

    emit:
    versions
    multiqc_files

}