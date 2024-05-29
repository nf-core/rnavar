//
// Subworkflow: Run GATK4 SplitNCigarReads with intervals, merge and index BAM file.
//

include { GATK4_SPLITNCIGARREADS } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { SAMTOOLS_MERGE         } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'

workflow SPLITNCIGAR {
    take:
    bam             // channel: [ val(meta), [ bam ], [bai] ]
    ch_fasta        // channel: [ fasta ]
    ch_fai          // channel: [ fai ]
    ch_dict         // channel: [ dict ]
    intervals       // channel: [ interval_list]

    main:
    ch_versions       = Channel.empty()

    bam_interval = bam.combine(intervals).map{ meta, bam, bai, intervals ->
        def new_meta = meta.clone()
        new_meta.id = meta.id + "_" + intervals.baseName
        new_meta.sample = meta.id
        [new_meta, bam, bai, intervals]
    }

    GATK4_SPLITNCIGARREADS(
        bam_interval,
        ch_fasta.map{ fasta -> [[id:'genome'], fasta] },
        ch_fai,
        ch_dict
    )
    bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

    bam_splitncigar_interval = bam_splitncigar
        .map{ meta, bam ->
            def new_meta = meta.clone()
            new_meta.id = meta.sample
            [new_meta, bam]
    }.groupTuple()

    SAMTOOLS_MERGE(bam_splitncigar_interval,
        ch_fasta,
        ch_fai.map{ fai -> [[id:fai.baseName], fai] })

    splitncigar_bam = SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX(splitncigar_bam)

    splitncigar_bam_bai = splitncigar_bam
        .join(SAMTOOLS_INDEX.out.bai, remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, remainder: true)
        .map{meta, bam, bai, csi ->
            if (bai) [meta, bam, bai]
            else [meta, bam, csi]
        }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        bam_bai     = splitncigar_bam_bai
        versions    = ch_versions
}
