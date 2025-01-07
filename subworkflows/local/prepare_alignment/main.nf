//
// Prepare input alignment files
//

include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    ch_cram     // [ val(meta), path(cram), path(crai) ]
    ch_bam      // [ val(meta), path(bam), path(bai) ]
    ch_fasta    // [ val(meta), path(fasta) ]
    ch_fai      // [ val(meta), path(fai) ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_CONVERT(
        ch_cram,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    def ch_bam_branch = ch_bam.branch { meta, bam, bai ->
        indexed: bai
            return [ meta, bam, bai ]
        not_indexed: !bai
            return [ meta, bam ]
    }

    def ch_bam_no_index = ch_bam_branch.not_indexed.mix(SAMTOOLS_CONVERT.out.cram)

    SAMTOOLS_INDEX(
        ch_bam_no_index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_bam_out = ch_bam_no_index
        .join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true)
        .mix(ch_bam_branch.indexed)

    emit:
    bam      = ch_bam_out // [ val(meta), path(bam), path(bai) ]
    versions = ch_versions
}
