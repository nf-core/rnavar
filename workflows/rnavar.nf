/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnavar.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.gtf,
    params.gff,
    params.dbsnp,
    params.dbsnp_tbi,
    params.known_indels,
    params.known_indels_tbi,
    params.snpeff_cache,
    params.vep_cache,
    params.star_index]

for (param in checkPathParamList) {if (param) file(param, checkIfExists: true)}

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if(!params.star_index && !params.gtf && !params.gff){ exit 1, "GTF|GFF3 file is required to build a STAR reference index! Use option -gtf|-gff to provide a GTF|GFF file." }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK                   } from '../subworkflows/local/input_check'              // Validate the input samplesheet.csv and prepare input channels
include { PREPARE_GENOME                } from '../subworkflows/local/prepare_genome'           // Build the genome index and other reference files
include { SPLITNCIGAR                   } from '../subworkflows/local/splitncigar'              // Splits reads that contain Ns in their cigar string
include { ANNOTATE                      } from '../subworkflows/local/annotate'                 // Annotate variants using snpEff or VEP or both

/*
========================================================================================
    IMPORT NF-CORE MODULES
========================================================================================
*/

include { FASTQC                        } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                       } from '../modules/nf-core/modules/multiqc/main'
include { CAT_FASTQ                     } from '../modules/nf-core/modules/cat/fastq/main'
include { GATK4_BASERECALIBRATOR        } from '../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_BEDTOINTERVALLIST       } from '../modules/nf-core/modules/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS       } from '../modules/nf-core/modules/gatk4/intervallisttools/main'
include { GATK4_HAPLOTYPECALLER         } from '../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS               } from '../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_INDEXFEATUREFILE        } from '../modules/nf-core/modules/gatk4/indexfeaturefile/main'
include { GATK4_VARIANTFILTRATION       } from '../modules/nf-core/modules/gatk4/variantfiltration/main'
include { SAMTOOLS_INDEX                } from '../modules/nf-core/modules/samtools/index/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    IMPORT NF-CORE SUBWORKFLOWS
========================================================================================
*/

include { ALIGN_STAR                    } from '../subworkflows/nf-core/align_star'         // Align reads to genome and sort and index the alignment file
include { MARKDUPLICATES                } from '../subworkflows/nf-core/markduplicates'     // Mark duplicates in the BAM file
include { RECALIBRATE                   } from '../subworkflows/nf-core/recalibrate'        // Estimate and correct systematic bias

/*
========================================================================================
    VARIABLES
========================================================================================
*/

// Check STAR alignment parameters
def prepareToolIndices  = params.aligner
def seq_platform        = params.seq_platform ? params.seq_platform : []
def seq_center          = params.seq_center ? params.seq_center : []

// Initialize file channels based on params
dbsnp                   = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()       : []
dbsnp_tbi               = params.dbsnp_tbi         ? Channel.fromPath(params.dbsnp_tbi).collect()   : []

// Initialize varaint annotation associated channels
def snpeff_db           = params.snpeff_db         ?: Channel.empty()
def vep_cache_version   = params.vep_cache_version ?: Channel.empty()
def vep_genome          = params.vep_genome        ?: Channel.empty()
def vep_species         = params.vep_species       ?: Channel.empty()
def snpeff_cache        = params.snpeff_cache      ? params.snpeff_cache : []
def vep_cache           = params.vep_cache         ? params.vep_cache : []

// MultiQC reporting
def multiqc_report = []

/*
========================================================================================
    RUN MAIN WORKFLOW RNAVAR
========================================================================================
*/

workflow RNAVAR {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (prepareToolIndices)
    ch_genome_bed = Channel.from([id:'genome.bed']).combine(PREPARE_GENOME.out.gene_bed)
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..meta.id.split('_').size()-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // PREPARE THE INTERVAL LIST FROM GTF FILE
    //

    ch_interval_list = Channel.empty()
    GATK4_BEDTOINTERVALLIST(ch_genome_bed, PREPARE_GENOME.out.dict)
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first().ifEmpty(null))

    //
    // MODULE: IntervalListTools from GATK4
    //

    ch_interval_list_split = Channel.empty()
    if (!params.skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(ch_interval_list)
        ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ meta, bed -> [bed] }.flatten()
    }
    else ch_interval_list_split = ch_interval_list

    //
    // SUBWORKFLOW: Alignment with STAR
    //

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()

    if (params.aligner == 'star') {
        ALIGN_STAR (
            ch_cat_fastq,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            seq_platform,
            seq_center
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        if (params.bam_csi_index) ch_genome_bam_index = ALIGN_STAR.out.csi
        ch_versions          = ch_versions.mix(ALIGN_STAR.out.versions.first().ifEmpty(null))

        // SUBWORKFLOW: Mark duplicates with Picard
        MARKDUPLICATES(ch_genome_bam)
        ch_genome_bam             = MARKDUPLICATES.out.bam.join(MARKDUPLICATES.out.bai, by: [0])
        ch_samtools_stats         = MARKDUPLICATES.out.stats
        ch_samtools_flagstat      = MARKDUPLICATES.out.flagstat
        ch_samtools_idxstats      = MARKDUPLICATES.out.idxstats
        ch_markduplicates_multiqc = MARKDUPLICATES.out.metrics
        if (params.bam_csi_index) ch_genome_bam_index = MARKDUPLICATES.out.csi
        ch_versions               = ch_versions.mix(MARKDUPLICATES.out.versions.first().ifEmpty(null))

        // Subworkflow - SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
        splitncigar_bam_bai = Channel.empty()
        SPLITNCIGAR(ch_genome_bam, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fai, PREPARE_GENOME.out.dict, ch_interval_list_split)
        splitncigar_bam_bai = SPLITNCIGAR.out.bam_bai
        ch_versions         = ch_versions.mix(SPLITNCIGAR.out.versions.first().ifEmpty(null))

        bam_variant_calling = Channel.empty()

        if(!params.skip_baserecalibration) {
            // MODULE: BaseRecalibrator from GATK4
            ch_bqsr_table = Channel.empty()
            known_sites     = Channel.from([params.dbsnp, params.known_indels]).collect()
            known_sites_tbi = Channel.from([params.dbsnp_tbi, params.known_indels_tbi]).collect()

            ch_interval_list_recalib = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
            splitncigar_bam_bai.combine(ch_interval_list_recalib)
                .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval]
            }.set{splitncigar_bam_bai_interval}

            GATK4_BASERECALIBRATOR(
                splitncigar_bam_bai_interval,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict,
                known_sites,
                known_sites_tbi
            )
            ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table
            ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first().ifEmpty(null))

            // MODULE: ApplyBaseRecalibrator from GATK4
            bam_applybqsr       = splitncigar_bam_bai.join(ch_bqsr_table, by: [0])
            bam_recalibrated_qc = Channel.empty()

            ch_interval_list_applybqsr = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
            bam_applybqsr.combine(ch_interval_list_applybqsr)
                .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]
            }.set{applybqsr_bam_bai_interval}

            RECALIBRATE(
                ('bamqc' in params.skip_qc),
                ('samtools' in params.skip_qc),
                applybqsr_bam_bai_interval,
                PREPARE_GENOME.out.dict,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.fasta
            )

            bam_variant_calling = RECALIBRATE.out.bam
            bam_recalibrated_qc = RECALIBRATE.out.qc
            ch_versions         = ch_versions.mix(RECALIBRATE.out.versions.first().ifEmpty(null))
        } else {
            bam_variant_calling = splitncigar_bam_bai
        }

        // MODULE: HaplotypeCaller from GATK4
        interval_flag = params.no_intervals
        haplotypecaller_vcf = Channel.empty()

        haplotypecaller_interval_bam = bam_variant_calling.combine(ch_interval_list_split)
            .map{ meta, bam, bai, interval_list ->
                new_meta = meta.clone()
                new_meta.id = meta.id + "_" + interval_list.baseName
                new_meta.sample = meta.id
                [new_meta, bam, bai, interval_list]
            }

        GATK4_HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict,
            dbsnp,
            dbsnp_tbi
        )

        haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
            .map{ meta, vcf ->
                meta.id = meta.sample
                [meta, vcf]}
            .groupTuple()

        ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first().ifEmpty(null))

        GATK4_MERGEVCFS(
            haplotypecaller_raw,
            PREPARE_GENOME.out.dict
        )
        haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
        ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first().ifEmpty(null))

        GATK4_INDEXFEATUREFILE(
            haplotypecaller_vcf
        )
        haplotypecaller_vcf     = haplotypecaller_vcf   //.join(GATK4_INDEXFEATUREFILE.out.index, by: [0])
        haplotypecaller_vcf_tbi = haplotypecaller_vcf.join(GATK4_INDEXFEATUREFILE.out.index, by: [0])
        ch_versions             = ch_versions.mix(GATK4_INDEXFEATUREFILE.out.versions.first().ifEmpty(null))
        final_vcf               = haplotypecaller_vcf

        // MODULE: VariantFiltration from GATK4
        if (!params.skip_variantfiltration) {

            GATK4_VARIANTFILTRATION(
                haplotypecaller_vcf_tbi,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict
            )

            filtered_vcf    = GATK4_VARIANTFILTRATION.out.vcf  //.join(GATK4_VARIANTFILTRATION.out.tbi, by: [0])
            final_vcf       = filtered_vcf
            ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first().ifEmpty(null))
        }

        if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff') || params.annotate_tools.contains('vep'))) {
            ANNOTATE(
                final_vcf,
                params.annotate_tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache)
            ch_versions = ch_versions.mix(ANNOTATE.out.versions.first().ifEmpty(null))
        }

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRnavar.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
