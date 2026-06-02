//
// Prepare input alignment files
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    bam // [ val(meta), path(bam), path(bai) ]
    cram // [ val(meta), path(cram), path(crai) ]

    main:
    def alignment_branch = bam
        .mix(cram)
        .branch { meta, reads, index ->
            indexed: index
            return [meta, reads, index]
            not_indexed: !index
            return [meta, reads]
        }

    SAMTOOLS_INDEX(alignment_branch.not_indexed)

    def alignment_out = alignment_branch.indexed
        .mix(alignment_branch.not_indexed)
        .join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)

    emit:
    reads_index = alignment_out // [ val(meta), path(bam|cram), path(bai|crai) ]
}
