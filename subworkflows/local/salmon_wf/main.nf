//
// RUN_SALMON
//

include { SALMON_INDEX  } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT  } from '../../../modules/nf-core/salmon/quant/main'


workflow SALMON_WF {
    take:
    samplesheet         // channel: [val(meta), reads ]
    fasta               // channel: [val(meta), fasta]
    salmon_index        // channel: [val(meta), index]
    transcriptome       // channel: [val(meta), fasta]
    gencode_gtf         // channel: [val(meta), gtf]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    if (!params.salmon_index){
        SALMON_INDEX (
            fasta.map{_meta,file -> file},
            transcriptome.map{_meta,file -> file}
        )
        salmon_index = SALMON_INDEX.out.index
    }

    SALMON_QUANT(
        samplesheet,
        salmon_index,
        gencode_gtf.map{_meta,file -> file},
        transcriptome.map{_meta,file -> file},
        false,
        'A'
    )
    multiqc_files = multiqc_files.mix(SALMON_QUANT.out.json_info)

    emit:
    versions
    multiqc_files

}