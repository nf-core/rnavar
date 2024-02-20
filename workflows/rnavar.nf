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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    params.star_index,
]

// only check if we are using the tools
if (params.annotate_tools && (params.annotate_tools.split(',').contains("snpeff") || params.annotate_tools.split(',').contains("merge"))) checkPathParamList.add(params.snpeff_cache)
if (params.annotate_tools && (params.annotate_tools.split(',').contains("vep")    || params.annotate_tools.split(',').contains("merge"))) checkPathParamList.add(params.vep_cache)

// Validate input parameters
WorkflowRnavar.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
    error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
}

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
include { BAM_MARKDUPLICATES } from '../subworkflows/local/bam_markduplicates/main'  // Mark duplicates in the BAM file
include { PREPARE_CACHE      } from '../subworkflows/local/prepare_cache/main'       // Download annotation cache if needed
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome/main'      // Build the genome index and other reference files
include { RECALIBRATE        } from '../subworkflows/local/recalibrate/main'         // Estimate and correct systematic bias
include { SPLITNCIGAR        } from '../subworkflows/local/splitncigar/main'         // Splits reads that contain Ns in their cigar string
include { VCF_ANNOTATE_ALL   } from '../subworkflows/local/vcf_annotate_all/main'    // Annotate variants using snpEff or VEP or both

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

// Initialize variant annotation associated channels
snpeff_db          = params.snpeff_db          ?: Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
if (params.snpeff_cache && params.annotate_tools && (params.annotate_tools.split(',').contains("snpeff") || params.annotate_tools.split(',').contains("merge"))) {
    def snpeff_annotation_cache_key = ''
    if (params.snpeff_cache == "s3://annotation-cache/snpeff_cache/") {
        snpeff_annotation_cache_key = "${params.snpeff_genome}.${params.snpeff_db}/"
    } else {
        snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
    }
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
    def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        if (params.snpeff_cache == "s3://annotation-cache/snpeff_cache/") {
            error("This path is not available within annotation-cache. Please check https://annotation-cache.github.io/ to create a request for it.")
        } else {
            error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
        }
    }
    snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
} else if (params.annotate_tools && (params.annotate_tools.split(',').contains("snpeff") || params.annotate_tools.split(',').contains("merge")) && !params.download_cache) {
        error("No cache for SnpEff or automatic download of said cache has been detected.\nPlease refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")
} else snpeff_cache = []

if (params.vep_cache && params.annotate_tools && (params.annotate_tools.split(',').contains("vep") || params.annotate_tools.split(',').contains("merge"))) {
    def vep_annotation_cache_key = ''
    if (params.vep_cache == "s3://annotation-cache/vep_cache/") {
        vep_annotation_cache_key = "${params.vep_cache_version}_${params.vep_genome}/"
    } else {
        vep_annotation_cache_key = params.use_annotation_cache_keys ? "${params.vep_cache_version}_${params.vep_genome}/" : ""
    }
    def vep_cache_dir = "${vep_annotation_cache_key}/${params.vep_species}"
    def vep_cache_path_full = file("$params.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        if (params.vep_cache == "s3://annotation-cache/vep_cache/") {
            error("This path is not available within annotation-cache. Please check https://annotation-cache.github.io/ to create a request for it.")
        } else {
            error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${params.vep_cache}.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
        }
    }
    vep_cache = Channel.fromPath(file("${params.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
} else if (params.annotate_tools && (params.annotate_tools.split(',').contains("vep") || params.annotate_tools.split(',').contains("merge")) && !params.download_cache) {
        error("No cache for VEP or automatic download of said cache has been detected.\nPlease refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")
} else vep_cache = []

vep_extra_files = []

if (params.dbnsfp && params.dbnsfp_tbi) {
    vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
    vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
}

if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
    vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
    vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
}

// MultiQC reporting
// def multiqc_report = []

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ch_exon_bed  = params.exon_bed ? Channel.fromPath(params.exon_bed)                                                       : Channel.empty()
ch_fasta_raw = params.fasta    ? Channel.fromPath(params.fasta).map{ fasta -> [ [ id:fasta.baseName ], fasta ] }.first() : Channel.empty()
ch_gff       = params.gff      ? Channel.fromPath(params.gff).first()                                                    : Channel.empty()
ch_gtf_raw   = params.gtf      ? Channel.fromPath(params.gtf).map{ gtf -> [ [ id:gtf.baseName ], gtf ] }.first()         : Channel.empty()

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

    ch_input = Channel.empty()

    if (params.input) ch_input = Channel.fromSamplesheet("input")

    ch_input = ch_input.map{ meta, fastq_1, fastq_2 ->
        if (fastq_2) return [ meta + [id: meta.sample], [ fastq_1, fastq_2 ] ]
        else return [ meta + [id: meta.sample], [ fastq_1 ] ]
    }.groupTuple()
    .branch { meta, fastq ->
        single  : fastq.size() == 1
            return [ meta + [single_end:true], fastq.flatten() ]
        multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
    }

    // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    snpeff_info     = params.snpeff_cache ? [] : Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])

    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info, snpeff_info)
        snpeff_cache = PREPARE_CACHE.out.snpeff_cache
        vep_cache    = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        ch_versions = ch_versions.mix(PREPARE_CACHE.out.versions)
    }

    // Prepare reference genome files

    PREPARE_GENOME(
        ch_fasta_raw,
        ch_gff,
        ch_gtf_raw,
        params.feature_type
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_genome_bed = params.exon_bed              ? Channel.fromPath(params.exon_bed).map{ it -> [ [id:'exon_bed'], it ] }.collect()
                                    : PREPARE_GENOME.out.exon_bed
    ch_dict       = params.dict                  ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : PREPARE_GENOME.out.dict

    ch_fasta      = params.fasta.endsWith('.gz') ? PREPARE_GENOME.out.fasta
                                    : ch_fasta_raw

    ch_fasta_fai  = params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()
                                    : PREPARE_GENOME.out.fasta_fai

    ch_gtf        = params.gtf.endsWith('.gz')   ? PREPARE_GENOME.out.gtf
                                    : ch_gtf_raw

    // MODULE: Concatenate FastQ files from same sample if required

    CAT_FASTQ(ch_input.multiple)

    ch_fastq = CAT_FASTQ.out.reads.mix(ch_input.single)

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // MODULE: Generate QC summary using FastQC
    FASTQC(ch_fastq)
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
            ch_fastq,
            PREPARE_GENOME.out.star_index.first(),
            PREPARE_GENOME.out.gtf.first(),
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
            ch_fasta,
            ch_fasta_fai.map{ it -> [[id:it.baseName], it] },
            ch_dict,
            ch_dbsnp.map{ it -> [[id:it.baseName], it] },
            ch_dbsnp_tbi.map{ it -> [[id:it.baseName], it] }
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
                    def new_meta = meta.clone()
                    new_meta.id = meta.sample
                    [new_meta, vcf]
                }.groupTuple()

            ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLERGVCF.out.versions.first().ifEmpty(null))
            //
            // MODULE: IndexFeatureFile from GATK4
            // Index the gVCF files
            //
            GATK4_INDEXFEATUREFILE(GATK4_HAPLOTYPECALLERGVCF.out.vcf)

            ch_haplotypecallergvcf_raw_index = GATK4_INDEXFEATUREFILE.out.index
                .map{ meta, idx ->
                    def new_meta = meta.clone()
                    new_meta.id = meta.sample
                    [new_meta, idx]
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

        //
        // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
        //
        if ((!params.skip_variantannotation) &&(params.annotate_tools) && (params.annotate_tools.split(',').contains('merge') || params.annotate_tools.split(',').contains('snpeff') || params.annotate_tools.split(',').contains('vep'))) {

            vep_fasta = (params.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]

            VCF_ANNOTATE_ALL(
                ch_final_vcf.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                vep_fasta,
                params.annotate_tools,
                params.snpeff_genome ? "${params.snpeff_genome}.${params.snpeff_db}" : "${params.genome}.${params.snpeff_db}",
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache,
                vep_extra_files)

            // Gather QC reports
            ch_reports  = ch_reports.mix(VCF_ANNOTATE_ALL.out.reports)
            ch_versions = ch_versions.mix(VCF_ANNOTATE_ALL.out.versions.first().ifEmpty(null))
        }

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
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
