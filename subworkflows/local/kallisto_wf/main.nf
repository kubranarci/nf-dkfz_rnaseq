//
// RUN_KALLISTO
//

include { KALLISTO_INDEX   } from '../../../modules/nf-core/kallisto/index/main'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant/main'
include { KALLISTO_RESCALE } from '../../../modules/local/kallisto_rescale/main'


workflow KALLISTO_WF {
    take:
    samplesheet         // channel: [val(meta), reads]
    index               // channel: [val(meta), star_index] 
    fasta               // channel: [val(meta), fasta]
    gencode_exclude_ch  // channel: [val(meta), gtf]
    kallisto_index      // channel: [val(meta), index]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    if (!params.kallisto_index) {
        KALLISTO_INDEX (
            fasta
        )

        kallisto_index = KALLISTO_INDEX.out.index
    }
    
    ch_kallisto_modes = channel.fromList( 
        params.run_kallisto ? params.run_kallisto.tokenize(',') : [] 
    )

    ch_kallisto_scatter = samplesheet
            .combine(ch_kallisto_modes)
            .map { meta, reads, mode ->
                def clean_mode = mode.trim()
                def new_meta = meta + [
                    id: "${meta.id}_${clean_mode}", 
                    strandedness: clean_mode
                ]
                return [ new_meta, reads ]
            }

    KALLISTO_QUANT(
        ch_kallisto_scatter,
        kallisto_index,
        [],[],"",""
    )
    multiqc_files = multiqc_files.mix(KALLISTO_QUANT.out.log.map{_meta, file -> file})

    KALLISTO_RESCALE(
       KALLISTO_QUANT.out.results,
       gencode_exclude_ch
    )

    emit:
    versions
    multiqc_files

}