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

    bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions)

    bam_splitncigar_interval = bam_splitncigar.map{ meta, bam_ -> [ meta + [id:meta.sample] - meta.subMap('sample'), bam_ ] }.groupTuple()

    SAMTOOLS_MERGE(bam_splitncigar_interval,
        fasta,
        fai.map{ fai_ -> [[id:fai_.baseName], fai_] })

    splitncigar_bam = SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    SAMTOOLS_INDEX(splitncigar_bam)

    splitncigar_bam_bai = splitncigar_bam
        .join(SAMTOOLS_INDEX.out.bai, remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, remainder: true) // TODO fix this bottleneck
        .map{meta, bam_, bai, csi ->
            if (bai) [meta, bam_, bai]
            else [meta, bam_, csi]
        }

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam_bai     = splitncigar_bam_bai
    versions    = ch_versions
}
