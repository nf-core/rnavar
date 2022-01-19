//
// PREPROCESSING OF ALIGNMENT FILE OVER THE INTERVALS IN PARALLEL
//

//params.gatk_splitncigar_options     = [:]
//params.samtools_index_options       = [:]
//params.samtools_merge_options       = [:]

include { GATK4_SPLITNCIGAR }   from '../nf-core/splitn_cigar_reads'                     //addParams(gatk_splitncigar_options: params.gatk_splitncigar_options, samtools_index_options: params.samtools_index_options)
include { SAMTOOLS_MERGE }      from '../../modules/nf-core/modules/samtools/merge/main' //addParams(options: params.samtools_merge_options)
include { SAMTOOLS_INDEX }      from '../../modules/nf-core/modules/samtools/index/main' //addParams(options: params.samtools_index_options)

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

    GATK4_SPLITNCIGAR(bam_interval, fasta, fasta_fai, fasta_dict)
    bam_splitncigar = GATK4_SPLITNCIGAR.out.bam
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGAR.out.versions.first())

    bam_splitncigar
        .map{ meta, bam ->
            meta.id = meta.sample
            [meta, bam]
    }.groupTuple().set{bam_splitncigar_interval}

    SAMTOOLS_MERGE(bam_splitncigar_interval, fasta)
    splitncigar_bam = SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX(splitncigar_bam)
    splitncigar_bam_bai = splitncigar_bam.join(SAMTOOLS_INDEX.out.bai)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        bam_bai     = splitncigar_bam_bai
        versions    = ch_versions
}
