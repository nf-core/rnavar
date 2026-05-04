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
    dict // channel: [mandatory] dict
    fai // channel: [mandatory] meta, fai
    fasta // channel: [mandatory] meta, fasta

    main:
    GATK4_APPLYBQSR(
        bam,
        fasta.map { _meta, fasta_ -> [fasta_] },
        fai.map { _meta, fai_ -> [fai_] },
        dict,
    )

    SAMTOOLS_INDEX(GATK4_APPLYBQSR.out.bam)

    def bam_recalibrated_index = GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)

    BAM_STATS_SAMTOOLS(bam_recalibrated_index, fasta.join(fai).collect())

    emit:
    bam = bam_recalibrated_index
}
