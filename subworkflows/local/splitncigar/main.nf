//
// Subworkflow: Run GATK4 SplitNCigarReads with intervals, merge and index BAM file.
//

include { GATK4_SPLITNCIGARREADS } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { SAMTOOLS_MERGE         } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'

workflow SPLITNCIGAR {
    take:
    bam          // channel: [ val(meta), [ bam ], [bai] ]
    fasta        // channel: [ val(meta), fasta ]
    fai          // channel: [ val(meta), fai ]
    dict         // channel: [ val(meta), dict ]
    intervals    // channel: [ interval_list]

    main:
    def ch_versions       = Channel.empty()

    def bam_interval = bam
        .combine(intervals)
        .map { meta, bam_, bai, intervals_ ->
            def new_meta = meta + [interval_count:intervals_ instanceof List ? intervals_.size() : 1]
            [ new_meta, bam_, bai, new_meta.interval_count > 1 ? intervals_ : [intervals_] ]
        }
        .transpose(by:3)
        .map { meta, bam_, bai, interval ->
            [ meta + [id:"${meta.id}_${interval.baseName}", sample: meta.id], bam_, bai, interval ]
        }

    GATK4_SPLITNCIGARREADS(bam_interval,
        fasta,
        fai,
        dict)

    def bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions)

    def bam_splitncigar_interval = bam_splitncigar
        .map{ meta, bam_ ->
            def new_meta = meta + [id:meta.sample] - meta.subMap('sample') - meta.subMap('interval_count')
            [ groupKey(new_meta, meta.interval_count), bam_ ]
        }
        .groupTuple()

    SAMTOOLS_MERGE(
        bam_splitncigar_interval,
        fasta,
        fai
    )

    def splitncigar_bam = SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    SAMTOOLS_INDEX(splitncigar_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    def splitncigar_bam_indices = SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.csi)
        .mix(SAMTOOLS_INDEX.out.crai)

    def splitncigar_bam_bai = splitncigar_bam
        .join(splitncigar_bam_indices, failOnDuplicate: true, failOnMismatch: true)

    emit:
    bam_bai     = splitncigar_bam_bai
    versions    = ch_versions
}
