/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

include { GATK4_APPLYBQSR as APPLYBQSR } from '../../modules/nf-core/modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX }               from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_STATS }               from '../../modules/nf-core/modules/samtools/stats/main'

workflow RECALIBRATE {
    take:
        skip_samtools  // boolean: true/false
        bam            // channel: [mandatory] bam
        dict           // channel: [mandatory] dict
        fai            // channel: [mandatory] fai
        fasta          // channel: [mandatory] fasta

    main:

    ch_versions = Channel.empty()

    bam_recalibrated_index = Channel.empty()
    bam_recalibrated       = Channel.empty()
    bam_reports            = Channel.empty()

    APPLYBQSR(bam, fasta, fai, dict)
    bam_recalibrated = APPLYBQSR.out.bam
    ch_versions = ch_versions.mix(APPLYBQSR.out.versions.first())

    SAMTOOLS_INDEX(bam_recalibrated)
    bam_recalibrated_index = bam_recalibrated
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map{meta, bam, bai, csi ->
            if (bai) [meta, bam, bai]
            else [meta, bam, csi]
        }

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    samtools_stats = Channel.empty()

    if (!skip_samtools) {
        SAMTOOLS_STATS(bam_recalibrated_index, [])
        samtools_stats = SAMTOOLS_STATS.out.stats
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    }
    bam_reports = samtools_stats


    emit:
        bam         = bam_recalibrated_index
        qc          = bam_reports

        versions    = ch_versions

}
