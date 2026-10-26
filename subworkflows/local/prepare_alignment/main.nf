//
// Prepare input alignment files
//

include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow PREPARE_ALIGNMENT {
    take:
    bam // [ val(meta), path(bam), path(bai) ]
    cram // [ val(meta), path(cram), path(crai) ]

    main:
    def input_reads = bam
        .mix(cram)
        .branch { meta, reads, index ->
            indexed: index
            return [meta, reads, index]
            not_indexed: !index && reads
            return [meta, reads]
        }

    SAMTOOLS_INDEX(input_reads.not_indexed)

    emit:
    reads_index = input_reads.indexed.mix(input_reads.not_indexed.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)) // [ val(meta), path(bam|cram), path(bai|crai) ]
}
