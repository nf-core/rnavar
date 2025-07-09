//
// Prepare input alignment files
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    cram // [ val(meta), path(cram), path(crai) ]
    bam  // [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()

    def alignment_branch = bam
        .mix(cram)
        .branch { meta, bam_, bai ->
            indexed: bai
            return [meta, bam_, bai]
            not_indexed_bam: !bai && bam_.extension == "bam"
            return [meta, bam_]
            not_indexed_cram: !bai && bam_.extension == "cram"
            return [meta, bam_]
        }

    def bam_no_index = alignment_branch.not_indexed_bam.mix(alignment_branch.not_indexed_cram)

    SAMTOOLS_INDEX(
        bam_no_index
    )

    def bam_indexed = alignment_branch.not_indexed_bam.join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true)
    def cram_indexed = alignment_branch.not_indexed_cram.join(SAMTOOLS_INDEX.out.crai, failOnMismatch: true, failOnDuplicate: true)

    def alignment_out = bam_indexed
        .mix(cram_indexed)
        .mix(alignment_branch.indexed)

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam      = alignment_out // [ val(meta), path(bam), path(bai) ]
    versions = ch_versions
}
