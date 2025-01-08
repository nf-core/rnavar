/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

include { GATK4_APPLYBQSR as APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS               } from '../../../modules/nf-core/samtools/stats/main'

workflow RECALIBRATE {
    take:
        skip_samtools  // boolean: true/false
        ch_bam         // channel: [mandatory] bam
        ch_dict        // channel: [mandatory] dict
        ch_fai         // channel: [mandatory] fai
        ch_fasta       // channel: [mandatory] fasta

    main:

    def ch_versions = Channel.empty()

    APPLYBQSR (
        ch_bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    def bam_recalibrated = APPLYBQSR.out.bam
    ch_versions = ch_versions.mix(APPLYBQSR.out.versions.first())

    SAMTOOLS_INDEX(bam_recalibrated)

    def bam_recalibrated_index = bam_recalibrated
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true) // TODO fix this bottleneck
        .map{meta, bam, bai, csi ->
            if (bai) [meta, bam, bai]
            else [meta, bam, csi]
        }

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def bam_reports = Channel.empty()

    if (!skip_samtools) {
        SAMTOOLS_STATS(bam_recalibrated_index, [[], []])
        bam_reports = SAMTOOLS_STATS.out.stats
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    }


    emit:
    bam         = bam_recalibrated_index
    qc          = bam_reports

    versions    = ch_versions

}
