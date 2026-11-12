/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

include { GATK4_APPLYBQSR    } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools'

workflow RECALIBRATE {
    take:
    bam // channel: [mandatory] bam
    fasta // channel: [mandatory] meta, fasta
    fai // channel: [mandatory] meta, fai
    dict // channel: [mandatory] meta, dict

    main:
    def fasta_fai_dict = fasta
        .combine(fai)
        .combine(dict)
        .map { meta_fasta, fasta_, _meta_fai, fai_, _meta_dict, dict_ -> [meta_fasta, fasta_, fai_, dict_]}
        .collect()

    GATK4_APPLYBQSR(bam, fasta_fai_dict, 'bam')

    SAMTOOLS_INDEX(GATK4_APPLYBQSR.out.bam)

    def bam_recalibrated_index = GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)

    BAM_STATS_SAMTOOLS(bam_recalibrated_index, fasta.combine(fai).map { meta, fasta_, _meta, fai_ -> [meta, fasta_, fai_] }.collect())

    emit:
    bam = bam_recalibrated_index
}
