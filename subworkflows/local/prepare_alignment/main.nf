//
// Prepare input alignment files
//

include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    ch_cram     // [ val(meta), path(cram), path(crai) ]
    ch_bam      // [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()

    def ch_alignment_branch = ch_bam
        .mix(ch_cram)
        .branch { meta, bam, bai ->
            indexed: bai
                return [ meta, bam, bai ]
            not_indexed_bam: !bai && bam.extension == "bam"
                return [ meta, bam ]
            not_indexed_cram: !bai && bam.extension == "cram"
                return [ meta, bam ]
        }

    def ch_bam_no_index = ch_alignment_branch.not_indexed_bam.mix(ch_alignment_branch.not_indexed_cram)

    SAMTOOLS_INDEX(
        ch_bam_no_index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_bam_indexed = ch_alignment_branch.not_indexed_bam
        .join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true)

    def ch_cram_indexed = ch_alignment_branch.not_indexed_cram
        .join(SAMTOOLS_INDEX.out.crai, failOnMismatch: true, failOnDuplicate: true)

    def ch_alignment_out = ch_bam_indexed
        .mix(ch_cram_indexed)
        .mix(ch_alignment_branch.indexed)

    emit:
    bam      = ch_alignment_out // [ val(meta), path(bam), path(bai) ]
    versions = ch_versions
}
