//
// RUN_RSEM
//

include { RSEM_PREPAREREFERENCE     } from '../../../modules/nf-core/rsem/preparereference/main'
include { RSEM_CALCULATEEXPRESSION  } from '../../../modules/nf-core/rsem/calculateexpression/main'

workflow RSEM_WF {
    take:
    transcipts          // channel [val(meta), bam]
    fasta               // channel: [val(meta), fasta]
    gencode_gtf_ch      // channel: [val(meta), gtf]
    rsem_index          // channel: [val(meta), index]
 
    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    if (!params.rsem_index){
        RSEM_PREPAREREFERENCE (
            fasta.map{_meta, file -> file},
            gencode_gtf_ch.map{_meta, file -> file}
        )
        rsem_index = RSEM_PREPAREREFERENCE.out.index

    }

    RSEM_CALCULATEEXPRESSION(
        transcipts,
        rsem_index
    )
    multiqc_files = multiqc_files.mix(RSEM_CALCULATEEXPRESSION.out.stat)

    emit:
    versions
    multiqc_files

}