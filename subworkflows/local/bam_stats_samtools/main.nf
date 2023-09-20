//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]

    main:
    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    SAMTOOLS_STATS(ch_bam_bai, [[],[]])

    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)
    ch_reports = ch_reports.mix(SAMTOOLS_FLAGSTAT.out.flagstat)
    ch_reports = ch_reports.mix(SAMTOOLS_IDXSTATS.out.idxstats)

    emit:
    reports = ch_reports
    versions = ch_versions                    // channel: [ versions.yml ]
}
