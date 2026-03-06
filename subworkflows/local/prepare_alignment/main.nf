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
            not_indexed_bam: !index && reads.extension == "bam"
            return [meta, reads]
            not_indexed_cram: !index && reads.extension == "cram"
            return [meta, reads]
        }

    SAMTOOLS_INDEX(alignment_branch.not_indexed_bam.mix(alignment_branch.not_indexed_cram))

    def alignment_out = alignment_branch.indexed
        .mix(alignment_branch.not_indexed_bam.join(SAMTOOLS_INDEX.out.bai, failOnMismatch: true, failOnDuplicate: true))
        .mix(alignment_branch.not_indexed_cram.join(SAMTOOLS_INDEX.out.crai, failOnMismatch: true, failOnDuplicate: true))

    emit:
    reads_index = alignment_out // [ val(meta), path(bam|cram), path(bai|crai) ]
}
