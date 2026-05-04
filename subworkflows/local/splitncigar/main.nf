//
// Subworkflow: Run GATK4 SplitNCigarReads with intervals, merge and index BAM file.
//

include { GATK4_SPLITNCIGARREADS } from '../../../modules/nf-core/gatk4/splitncigarreads'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_MERGE         } from '../../../modules/nf-core/samtools/merge'

workflow SPLITNCIGAR {
    take:
    bam // channel: [ val(meta), [ bam ], [bai] ]
    fasta // channel: [ val(meta), fasta ]
    fai // channel: [ val(meta), fai ]
    dict // channel: [ val(meta), dict ]
    intervals // channel: [ interval_list]

    main:
    def bam_for_splincigarreads = channel.empty()

    if (intervals) {
        bam_for_splincigarreads = bam
            .combine(intervals)
            .map { meta, bam_, bai, intervals_ ->
                [
                    meta + [interval_count: intervals_ instanceof List ? intervals_.size() : 1],
                    bam_,
                    bai,
                    intervals_ instanceof List ? intervals_ : [intervals_],
                ]
            }
            .transpose(by: 3)
            .map { meta, bam_, bai, interval -> [meta + [id: "${meta.id}_${interval.baseName}", sample: meta.id], bam_, bai, interval] }
    }
    else {
        bam_for_splincigarreads = bam.map { meta, bam_, bai -> [meta + [interval_count: 1, sample: meta.id], bam_, bai, []] }
    }

    GATK4_SPLITNCIGARREADS(bam_for_splincigarreads, fasta, fai, dict)

    def bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam

    def bam_to_merge = bam_splitncigar
        .map { meta, bam_ -> [groupKey(meta + [id: meta.sample] - meta.subMap('sample') - meta.subMap('interval_count'), meta.interval_count), bam_] }
        .groupTuple()
        .map { meta, bam_ -> [meta, bam_, []] }

    SAMTOOLS_MERGE(
        bam_to_merge,
        fasta.join(fai).map { meta, _fasta, _fai -> [meta, _fasta, _fai, []] }.collect(),
    )

    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)

    emit:
    bam_bai = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_INDEX.out.index, failOnDuplicate: true, failOnMismatch: true)
}
