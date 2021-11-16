//
// GATK4 SplitNCigarReads
//

params.gatk_splitncigar_options = [:]
params.samtools_index_options   = [:]
params.samtools_stats_options   = [:]

include { GATK4_SPLITNCIGARREADS } from '../../modules/nf-core/modules/gatk4/splitncigarreads/main' addParams( options: params.gatk_splitncigar_options )
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/modules/samtools/index/main'         addParams( options: params.samtools_index_options )

workflow GATK4_SPLITNCIGAR {
    take:
    bam             // channel: [ val(meta), [ bam ] ]
    fasta           // channel: [ fasta ]
    fasta_fai      // channel: [ fai ]
    fasta_dict     // channel: [ dict ]

    main:

    ch_versions = Channel.empty()

    //
    // GATK4
    //
    GATK4_SPLITNCIGARREADS ( bam, fasta, fasta_fai, fasta_dict)
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( GATK4_SPLITNCIGARREADS.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    GATK4_SPLITNCIGARREADS.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }


    emit:
    bam              = GATK4_SPLITNCIGARREADS.out.bam         // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    csi              = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), [ csi ] ]
    versions         = ch_versions                       // channel: [ versions.yml ]
}
