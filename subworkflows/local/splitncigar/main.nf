//
// Subworkflow: Run GATK4 SplitNCigarReads with intervals, merge and index BAM file.
//

include { GATK4_SPLITNCIGARREADS } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { SAMTOOLS_MERGE         } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'

workflow SPLITNCIGAR {
    take:
    bam          // channel: [ val(meta), [ bam ], [bai] ]
    fasta        // channel: [ fasta ]
    fai          // channel: [ fai ]
    dict         // channel: [ dict ]
    intervals    // channel: [ interval_list]

    main:
    def ch_versions       = Channel.empty()

    def bam_interval = bam.combine(intervals).map{ meta, bam_, bai, intervals_ -> [ meta + [sample:meta.id], bam_, bai, intervals_ ] }

    GATK4_SPLITNCIGARREADS(bam_interval,
        fasta,
        fai.map{ fai_ -> [[id:'genome'], fai_] },
        dict)

    def bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions)

    def bam_splitncigar_interval = bam_splitncigar.map{ meta, bam_ -> [ meta + [id:meta.sample] - meta.subMap('sample'), bam_ ] }.groupTuple()

    SAMTOOLS_MERGE(
        bam_splitncigar_interval,
        fasta,
        fai.map{ fai_ -> [[id:fai_.baseName], fai_] }
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
