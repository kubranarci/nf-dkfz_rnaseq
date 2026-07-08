//
// RUN_STAR_ARRIBA
//

include { STAR_ALIGN          } from '../../../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'
include { SAMTOOLS_SORT       } from '../../../modules/nf-core/samtools/sort/main'
include { SAMBAMBA_MARKDUP    } from '../../../modules/nf-core/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT   } from '../../../modules/nf-core/sambamba/flagstat/main'
include { ARRIBA_DOWNLOAD     } from '../../../modules/nf-core/arriba/download/main'
include { FIND_CONCATENATE    } from '../../../modules/nf-core/find/concatenate/main' 
include { SAMTOOLS_VIEW       } from '../../../modules/nf-core/samtools/view/main'
include { GUNZIP              } from '../../../modules/nf-core/gunzip/main'
workflow STAR_ARRIBA_WF {
    take:
    samplesheet     // channel: [val(meta), reads ]
    arriba          // channel: [val(meta), arriba_index] 
    fasta           // channel: [val(meta), fasta]
    fai             // channel: [val(meta), fasta]
    fasta_with_virus // channel: [val(meta), fasta]
    gtf             // channel: [val(meta), gtf]
    viral           // channel: [val(meta), fasta]
    blacklist       // channel: [file]
    known_fusions   // channel: [file]
    cytobands       // channel: [file]
    protein_domains // channel: [file]

    main:

    versions = channel.empty()
    multiqc_files = channel.empty()

    // Flag helper to ensure the string contains real data
    def has_virus_list = params.ignoreViruses != null && params.ignoreViruses != "" && params.ignoreViruses != "false"

    // 1. Download missing assets if a virus list is specified
    if ((!params.fasta_with_virus || !params.blacklist || !params.cytobands || !params.knownFusions || !params.proteinDomains) && has_virus_list) {
        ARRIBA_DOWNLOAD(
            params.genome
        )
        blacklist       = ARRIBA_DOWNLOAD.out.blacklist
        known_fusions   = ARRIBA_DOWNLOAD.out.known_fusions
        cytobands       = ARRIBA_DOWNLOAD.out.cytobands
        protein_domains = ARRIBA_DOWNLOAD.out.protein_domains
        viral           = ARRIBA_DOWNLOAD.out.viral_genome
    }



    // 2. Concatenate Human + Viral FASTAs if the combined file is missing
    if (!params.fasta_with_virus && has_virus_list) {
        GUNZIP(
            viral.map{file -> tuple ([id:'viral'], file)}
        )
        ch_fasta_to_concat = fasta
            .combine(GUNZIP.out.gunzip)
            .map { meta, human_fa, _meta2, viral_fa -> 
                tuple(meta, [ human_fa, viral_fa ]) 
            }

        FIND_CONCATENATE(
            ch_fasta_to_concat
        )
        fasta_with_virus = FIND_CONCATENATE.out.file_out
    }

    // 3. Generate the Arriba-optimized STAR index if missing
    if (!params.star_arriba_index && has_virus_list) {
        STAR_GENOMEGENERATE (
            fasta_with_virus,
            gtf
        )
        arriba = STAR_GENOMEGENERATE.out.index
    }

    // 4. Align Reads (Always runs)
    STAR_ALIGN (
        samplesheet,
        arriba,
        [[],[]],
        true,
        params.seq_platform,
        params.seq_center
    )
    versions = STAR_ALIGN.out.versions
    multiqc_files = multiqc_files.mix(STAR_ALIGN.out.log_final)

    ch_out_bam  = STAR_ALIGN.out.bam

    if (has_virus_list) {
        SAMTOOLS_VIEW(
            ch_out_bam.map{meta, file -> tuple(meta, file, [])},
            fasta.combine(fai).map{meta, fasta_ch, _meta2, fai_ch -> tuple(meta, fasta_ch, fai_ch)},
            [],
            "bai"
        )
        ch_out_bam = SAMTOOLS_VIEW.out.bam
    }

    emit:
    versions
    multiqc_files
    ch_out_bam
    arriba
    blacklist      
    known_fusions 
    cytobands     
    protein_domains
    fasta_with_virus 
}