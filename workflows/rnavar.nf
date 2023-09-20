/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

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
    params.star_index,
]

for(param in checkPathParamList) {if (param) file(param, checkIfExists: true)}

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (!params.star_index && !params.gtf && !params.gff){ exit 1, "GTF|GFF3 file is required to build a STAR reference index! Use option --gtf|--gff to provide a GTF|GFF file." }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GTF2BED            } from '../modules/local/gtf2bed/main'

include { ALIGN_STAR         } from '../subworkflows/local/align_star/main'          // Align reads to genome and sort and index the alignment file
include { ANNOTATE           } from '../subworkflows/local/annotate/main'            // Annotate variants using snpEff or VEP or both
include { BAM_MARKDUPLICATES } from '../subworkflows/local/bam_markduplicates/main'  // Mark duplicates in the BAM file
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome/main'      // Build the genome index and other reference files
include { RECALIBRATE        } from '../subworkflows/local/recalibrate/main'         // Estimate and correct systematic bias
include { SPLITNCIGAR        } from '../subworkflows/local/splitncigar/main'         // Splits reads that contain Ns in their cigar string

/*
========================================================================================
    IMPORT NF-CORE MODULES
========================================================================================
*/

include { FASTQC                                             } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                            } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                                          } from '../modules/nf-core/cat/fastq/main'
include { GATK4_BASERECALIBRATOR                             } from '../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_BEDTOINTERVALLIST                            } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS                            } from '../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_HAPLOTYPECALLER                              } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_HAPLOTYPECALLER as GATK4_HAPLOTYPECALLERGVCF } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS                                    } from '../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_COMBINEGVCFS                                 } from '../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_INDEXFEATUREFILE                             } from '../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_VARIANTFILTRATION                            } from '../modules/nf-core/gatk4/variantfiltration/main'
include { SAMTOOLS_INDEX                                     } from '../modules/nf-core/samtools/index/main'
include { TABIX_TABIX as TABIX                               } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIXGVCF                           } from '../modules/nf-core/tabix/tabix/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                        } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
========================================================================================
    VARIABLES
========================================================================================
*/

// // Check STAR alignment parameters
// def prepareToolIndices  = params.aligner
def seq_platform = params.seq_platform ?: []
def seq_center   = params.seq_center   ?: []

// // Initialize file channels based on params
ch_dbsnp                = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()               : Channel.empty()
ch_dbsnp_tbi            = params.dbsnp_tbi         ? Channel.fromPath(params.dbsnp_tbi).collect()           : Channel.empty()
ch_known_indels         = params.known_indels      ? Channel.fromPath(params.known_indels).collect()        : Channel.empty()
ch_known_indels_tbi     = params.known_indels_tbi  ? Channel.fromPath(params.known_indels_tbi).collect()    : Channel.empty()

// // Initialize variant annotation associated channels
// ch_snpeff_db            = params.snpeff_db         ?:   Channel.empty()
// ch_vep_cache_version    = params.vep_cache_version ?:   Channel.empty()
// ch_vep_genome           = params.vep_genome        ?:   Channel.empty()
// ch_vep_species          = params.vep_species       ?:   Channel.empty()
// ch_snpeff_cache         = params.snpeff_cache      ?    Channel.fromPath(params.snpeff_cache).collect()  : []
// ch_vep_cache            = params.vep_cache         ?    Channel.fromPath(params.vep_cache).collect()     : []

// MultiQC reporting
// def multiqc_report = []



// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ch_exon_bed = params.exon_bed ? Channel.fromPath(params.exon_bed)                                                       : Channel.empty()
ch_fasta    = params.fasta    ? Channel.fromPath(params.fasta).map{ fasta -> [ [ id:fasta.baseName ], fasta ] }.first() : Channel.empty()
ch_gff      = params.gff      ? Channel.fromPath(params.gff).first()                                                    : Channel.empty()
ch_gtf      = params.gtf      ? Channel.fromPath(params.gtf).map{ gtf -> [ [ id:gtf.baseName ], gtf ] }.first()         : Channel.empty()

/*
========================================================================================
    RUN MAIN WORKFLOW RNAVAR
========================================================================================
*/

workflow RNAVAR {
    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()

    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

    ch_from_samplesheet = Channel.empty()

    if (params.input) ch_from_samplesheet = Channel.fromSamplesheet("input")

    ch_fastq = ch_from_samplesheet.map{ meta, fastq_1, fastq_2 ->
        if (fastq_2) return [ meta + [id: meta.sample], [ fastq_1, fastq_2 ] ]
        else return [ meta + [id: meta.sample], [ fastq_1 ] ]
    }.groupTuple()
    .branch { meta, fastq ->
        single  : fastq.size() == 1
            return [ meta, fastq.flatten() ]
        multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
    }

    // Prepare reference genome files

    PREPARE_GENOME(
        ch_fasta,
        ch_gff,
        ch_gtf,
        params.feature_type
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_genome_bed = params.exon_bed  ? Channel.fromPath(params.exon_bed).map{ it -> [ [id:'exon_bed'], it ] }.collect()
                                    : PREPARE_GENOME.out.exon_bed
    ch_dict       = params.dict      ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : ch_dict
    ch_fasta_fai  = params.fasta_fai ? Channel.fromPath(params.fasta_fai).collect()
                                    : PREPARE_GENOME.out.fasta_fai

    // MODULE: Concatenate FastQ files from same sample if required

    CAT_FASTQ(ch_fastq.multiple)

    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // MODULE: Generate QC summary using FastQC
    FASTQC(ch_cat_fastq)
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
    //

    GATK4_BEDTOINTERVALLIST(ch_genome_bed, ch_dict)
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first().ifEmpty(null))

    //
    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    //
    ch_interval_list_split = Channel.empty()
    if (!params.skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(ch_interval_list)
        ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ meta, bed -> [bed] }.flatten()
    }
    else ch_interval_list_split = ch_interval_list

    //
    // SUBWORKFLOW: Perform read alignment using STAR aligner
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
        ALIGN_STAR(
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

        // Gather QC reports
        ch_reports           = ch_reports.mix(ALIGN_STAR.out.reports)
        ch_reports           = ch_reports.mix(ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_versions          = ch_versions.mix(ALIGN_STAR.out.versions.first().ifEmpty(null))

        //
        // SUBWORKFLOW: Mark duplicates with GATK4
        //
        BAM_MARKDUPLICATES(
            ch_genome_bam,
            ch_fasta.map{ meta, fasta -> [fasta] },
            ch_fasta_fai,
            [])

        ch_genome_bam             = BAM_MARKDUPLICATES.out.bam

        //Gather QC reports
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES.out.reports.collect{it[1]}.ifEmpty([]))
        ch_versions               = ch_versions.mix(BAM_MARKDUPLICATES.out.versions.first().ifEmpty(null))

        //
        // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string(e.g. spanning splicing events in RNAseq data).
        //
        ch_splitncigar_bam_bai = Channel.empty()
        SPLITNCIGAR(
            ch_genome_bam,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_interval_list_split
        )
        ch_splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
        ch_versions             = ch_versions.mix(SPLITNCIGAR.out.versions.first().ifEmpty(null))

        //
        // MODULE: BaseRecalibrator from GATK4
        // Generates a recalibration table based on various co-variates
        //
        ch_bam_variant_calling = Channel.empty()
        if (!params.skip_baserecalibration) {
            ch_bqsr_table   = Channel.empty()
            // known_sites is made by grouping both the dbsnp and the known indels ressources
            // they can either or both be optional
            ch_known_sites     = ch_dbsnp.concat(ch_known_indels).collect()
            ch_known_sites_tbi = ch_dbsnp_tbi.concat(ch_known_indels_tbi).collect()

            ch_interval_list_recalib = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
            ch_splitncigar_bam_bai_interval = ch_splitncigar_bam_bai.combine(ch_interval_list_recalib)
                .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval] }

            GATK4_BASERECALIBRATOR(
                ch_splitncigar_bam_bai_interval,
                ch_fasta.map{ meta, fasta -> [fasta] },
                ch_fasta_fai,
                ch_dict.map{ meta, dict -> [dict] },
                ch_known_sites,
                ch_known_sites_tbi
            )
            ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table

            // Gather QC reports
            ch_reports  = ch_reports.mix(ch_bqsr_table.map{ meta, table -> table})
            ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first().ifEmpty(null))

            ch_bam_applybqsr       = ch_splitncigar_bam_bai.join(ch_bqsr_table)
            ch_bam_recalibrated_qc = Channel.empty()

            ch_interval_list_applybqsr = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
            ch_bam_applybqsr.combine(ch_interval_list_applybqsr)
                .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]}
                .set{ch_applybqsr_bam_bai_interval}

            //
            // MODULE: ApplyBaseRecalibrator from GATK4
            // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
            //
            RECALIBRATE(
                params.skip_multiqc,
                ch_applybqsr_bam_bai_interval,
                ch_dict.map{ meta, dict -> [dict] },
                ch_fasta_fai,
                ch_fasta.map{ meta, fasta -> [fasta] }
            )

            ch_bam_variant_calling = RECALIBRATE.out.bam
            ch_bam_recalibrated_qc = RECALIBRATE.out.qc

            // Gather QC reports
            ch_reports  = ch_reports.mix(RECALIBRATE.out.qc.collect{it[1]}.ifEmpty([]))
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions.first().ifEmpty(null))
        } else {
            ch_bam_variant_calling = ch_splitncigar_bam_bai
        }

        interval_flag = params.no_intervals
        // Run haplotyper even in the absence of dbSNP files
        if (!params.dbsnp){
            ch_dbsnp = []
            ch_dbsnp_tbi = []
        }

        ch_haplotypecaller_vcf = Channel.empty()
        ch_haplotypecaller_interval_bam = ch_bam_variant_calling.combine(ch_interval_list_split)
            .map{ meta, bam, bai, interval_list ->
                [meta + [id:meta.id + "_" + interval_list.baseName], bam, bai, interval_list, []]
            }

        //
        // MODULE: HaplotypeCaller from GATK4
        // Calls germline SNPs and indels via local re-assembly of haplotypes.
        //

        GATK4_HAPLOTYPECALLER(
            ch_haplotypecaller_interval_bam,
            ch_fasta.map{ meta, fasta -> [fasta] },
            ch_fasta_fai,
            ch_dict.map{ meta, dict -> [dict] },
            ch_dbsnp,
            ch_dbsnp_tbi
        )


        ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
            .map{ meta, vcf ->
                meta.id = meta.sample
                [meta, vcf]}
            .groupTuple()

        ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first().ifEmpty(null))

        //
        // MODULE: MergeVCFS from GATK4
        // Merge multiple VCF files into one VCF
        //
        GATK4_MERGEVCFS(ch_haplotypecaller_raw, ch_dict)
        ch_haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
        ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first().ifEmpty(null))

        if (params.generate_gvcf){
            GATK4_HAPLOTYPECALLERGVCF(
                ch_haplotypecaller_interval_bam,
                ch_fasta.map{ meta, fasta -> [fasta] },
                ch_fasta_fai,
                ch_dict.map{ meta, dict -> [dict] },
                ch_dbsnp,
                ch_dbsnp_tbi
            )

            ch_haplotypecallergvcf_raw = GATK4_HAPLOTYPECALLERGVCF.out.vcf
                .map{ meta, vcf ->
                    meta.id = meta.sample
                    [meta, vcf]
                }.groupTuple()

            ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLERGVCF.out.versions.first().ifEmpty(null))
            //
            // MODULE: IndexFeatureFile from GATK4
            // Index the gVCF files
            //
            GATK4_INDEXFEATUREFILE(GATK4_HAPLOTYPECALLERGVCF.out.vcf)

            ch_haplotypecallergvcf_raw_index = GATK4_INDEXFEATUREFILE.out.index
                .map{ meta, idx ->
                    meta.id = meta.sample
                    [meta, idx]
                }.groupTuple()

            ch_versions  = ch_versions.mix(GATK4_INDEXFEATUREFILE.out.versions.first().ifEmpty(null))

            // MODULE: CombineGVCFS from GATK4
            // Merge multiple GVCF files into one GVCF

            ch_haplotypecallergvcf_raw_tbi = ch_haplotypecallergvcf_raw
                .join(ch_haplotypecallergvcf_raw_index, remainder: true)

            GATK4_COMBINEGVCFS(
                ch_haplotypecallergvcf_raw_tbi,
                ch_fasta.map{ meta, fasta -> [fasta] },
                ch_fasta_fai,
                ch_dict.map{ meta, dict -> [dict] }
            )
            ch_haplotypecaller_gvcf = GATK4_COMBINEGVCFS.out.combined_gvcf
            ch_versions  = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions.first().ifEmpty(null))

            //
            // MODULE: Index the VCF using TABIX
            //
            TABIXGVCF(ch_haplotypecaller_gvcf)

            ch_haplotypecaller_gvcf_tbi = ch_haplotypecaller_gvcf
                .join(TABIXGVCF.out.tbi, remainder: true)
                .join(TABIXGVCF.out.csi, remainder: true)
                .map{meta, vcf, tbi, csi ->
                    if (tbi) [meta, vcf, tbi]
                    else [meta, vcf, csi]
                }

            ch_versions  = ch_versions.mix(TABIXGVCF.out.versions.first().ifEmpty(null))

        }

        //
        // MODULE: Index the VCF using TABIX
        //
        TABIX(ch_haplotypecaller_vcf)

        ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_vcf
            .join(TABIX.out.tbi, remainder: true)
            .join(TABIX.out.csi, remainder: true)
            .map{meta, vcf, tbi, csi ->
                if (tbi) [meta, vcf, tbi]
                else [meta, vcf, csi]
            }

        ch_versions     = ch_versions.mix(TABIX.out.versions.first().ifEmpty(null))
        ch_final_vcf    = ch_haplotypecaller_vcf

        //
        // MODULE: VariantFiltration from GATK4
        // Filter variant calls based on certain criteria
        //
        if (!params.skip_variantfiltration && !params.bam_csi_index ) {

            GATK4_VARIANTFILTRATION(
                ch_haplotypecaller_vcf_tbi,
                ch_fasta,
                ch_fasta_fai.map{ it -> [ [id:'fai'], it ] },
                ch_dict
            )

            ch_filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
            ch_final_vcf    = ch_filtered_vcf
            ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first().ifEmpty(null))
        }

    //     //
    //     // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
    //     //
    //     if ((!params.skip_variantannotation) &&(params.annotate_tools) &&(params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff') || params.annotate_tools.contains('vep'))) {
    //         ANNOTATE(
    //             ch_final_vcf,
    //             params.annotate_tools,
    //             ch_snpeff_db,
    //             ch_snpeff_cache,
    //             ch_vep_genome,
    //             ch_vep_species,
    //             ch_vep_cache_version,
    //             ch_vep_cache)

    //         // Gather QC reports
    //         ch_reports  = ch_reports.mix(ANNOTATE.out.reports)
    //         ch_versions = ch_versions.mix(ANNOTATE.out.versions.first().ifEmpty(null))
    //     }

    }

    ch_version_yaml = Channel.empty()
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

    //
    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    //
    if (!params.skip_multiqc){
        workflow_summary    = WorkflowRnavar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRnavar.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_reports.view()

        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix(ch_version_yaml)
        multiqc_files = multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        multiqc_files = multiqc_files.mix(ch_reports.collect().ifEmpty([]))

        MULTIQC(multiqc_files.collect(), ch_multiqc_config.collect().ifEmpty([]), ch_multiqc_custom_config.collect().ifEmpty([]), ch_multiqc_logo.collect().ifEmpty([]))

        multiqc_report = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
