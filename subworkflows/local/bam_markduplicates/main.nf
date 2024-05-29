//
// MARKDUPLICATES AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools'
include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX       } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES {
    take:
    bam                    // channel: [mandatory] [ meta, bam ]
    fasta                  // channel: [mandatory] [ fasta ]
    fasta_fai              // channel: [mandatory] [ fasta_fai ]
    intervals_bed_combined // channel: [optional]  [ intervals_bed ]

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // RUN MARKUPDUPLICATES
    GATK4_MARKDUPLICATES(bam, fasta, fasta_fai)

    SAMTOOLS_INDEX(GATK4_MARKDUPLICATES.out.bam)

    ch_bam_index = GATK4_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai, remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, remainder: true)
        .map{meta, bam, bai, csi ->
            if (bai) [meta, bam, bai]
            else [meta, bam, csi]
        }

    BAM_STATS_SAMTOOLS(ch_bam_index, fasta)

    // Gather all reports generated
    ch_reports = ch_reports.mix(GATK4_MARKDUPLICATES.out.metrics)
    ch_reports = ch_reports.mix(BAM_STATS_SAMTOOLS.out.stats)
    ch_reports = ch_reports.mix(BAM_STATS_SAMTOOLS.out.flagstat)
    ch_reports = ch_reports.mix(BAM_STATS_SAMTOOLS.out.idxstats)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam     = ch_bam_index
    reports = ch_reports

    versions = ch_versions    // channel: [ versions.yml ]
}
