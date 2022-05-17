//
// Subworkflow: Run GATK4 SplitNCigarReads with intervals, merge and index BAM file.
//

include { GATK4_SPLITNCIGARREADS } from '../../modules/nf-core/modules/gatk4/splitncigarreads/main'
include { SAMTOOLS_MERGE         } from '../../modules/nf-core/modules/samtools/merge/main'
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/modules/samtools/index/main'

workflow SPLITNCIGAR {
    take:
    bam             // channel: [ val(meta), [ bam ], [bai] ]
    fasta           // channel: [ fasta ]
    fasta_fai       // channel: [ fai ]
    fasta_dict      // channel: [ dict ]
    intervals       // channel: [ interval_list]

    main:

    ch_versions       = Channel.empty()

    bam.combine(intervals)
        .map{ meta, bam, bai, intervals ->
        new_meta = meta.clone()
        new_meta.id = meta.id + "_" + intervals.baseName
        new_meta.sample = meta.id
        [new_meta, bam, bai, intervals]
    }.set{bam_interval}

    GATK4_SPLITNCIGARREADS(bam_interval, fasta, fasta_fai, fasta_dict)
    bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

    bam_splitncigar
        .map{ meta, bam ->
            new_meta = meta.clone()
            new_meta.id = meta.sample
            [new_meta, bam]
    }.groupTuple().set{bam_splitncigar_interval}

    SAMTOOLS_MERGE(bam_splitncigar_interval, fasta)
    splitncigar_bam = SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX(splitncigar_bam)
    splitncigar_bam_bai = splitncigar_bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map{meta, bam, bai, csi ->
            if (bai) [meta, bam, bai]
            else [meta, bam, csi]
        }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        bam_bai     = splitncigar_bam_bai
        versions    = ch_versions
}
