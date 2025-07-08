/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_INDEX  } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_STATS  } from '../../../modules/nf-core/samtools/stats'

workflow RECALIBRATE {
    take:
    skip_samtools // boolean: true/false
    bam           // channel: [mandatory] bam
    dict          // channel: [mandatory] dict
    fai           // channel: [mandatory] fai
    fasta         // channel: [mandatory] fasta

    main:

    def ch_reports = Channel.empty()
    def ch_versions = Channel.empty()

    GATK4_APPLYBQSR(
        bam,
        fasta,
        fai,
        dict,
    )

    SAMTOOLS_INDEX(GATK4_APPLYBQSR.out.bam)

    def bam_indices = SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.csi)
        .mix(SAMTOOLS_INDEX.out.crai)

    def bam_recalibrated_index = GATK4_APPLYBQSR.out.bam.join(bam_indices, failOnMismatch: true, failOnDuplicate: true)

    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    if (!skip_samtools) {
        SAMTOOLS_STATS(bam_recalibrated_index, [[], []])
        ch_reports = SAMTOOLS_STATS.out.stats
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    }

    emit:
    bam      = bam_recalibrated_index
    qc       = ch_reports
    versions = ch_versions
}
