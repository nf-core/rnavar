/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnavar.initialise(params, log)

// Check input path parameters to see if they exist
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.gtf,
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

if(!params.star_index and !params.gtf){ exit 1, "GTF file is required to build a STAR reference index! Use option -gtf to provide a GTF file." }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'  addParams( options: [publish_files : ['_versions.yml':'']] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Check alignment parameters
def prepareToolIndices  = []
prepareToolIndices = params.aligner

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]
def untar_options          = [publish_files: false]
def samtools_sort_genome_options = modules['samtools_sort_genome']

if (!params.save_reference) modules['star_genomegenerate']['publish_files'] = false

modules['samtools_index_genome'].args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

if (params.save_align_intermeds) {
    samtools_sort_genome_options.publish_files.put('bam','')
    samtools_sort_genome_options.publish_files.put('bai','')
    samtools_sort_genome_options.publish_files.put('csi','')
}

if (!params.save_merged_fastq) modules['cat_fastq']['publish_files'] = false

// Additional paramaters to modules based on params
modules['star_align'].args += params.save_unaligned ? Utils.joinModuleArgs(['--outReadsUnmapped Fastx']) : ''
modules['star_align'].args += params.star_twopass ? Utils.joinModuleArgs(['--twopassMode Basic']) : ''
if (params.save_align_intermeds) modules['star_align'].publish_files.put('bam','')
if (params.save_unaligned)       modules['star_align'].publish_files.put('fastq.gz','unmapped')

modules['picard_markduplicates_samtools'].args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''
modules['picard_markduplicates'].args += params.remove_duplicates ? Utils.joinModuleArgs(['REMOVE_DUPLICATES=true']) : ''
modules['gatk_intervallisttools'].args += Utils.joinModuleArgs(["--SCATTER_COUNT $params.scatter_count"])

modules['gatk_haplotypecaller'].args += Utils.joinModuleArgs(["--standard-min-confidence-threshold-for-calling $params.stand_call_conf"])

if (params.window)    modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--window $params.window"])
if (params.cluster)   modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--cluster $params.cluster"])
if (params.fs_filter) modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--filter-name \"FS\" --filter \"FS > $params.fs_filter\" "])
if (params.qd_filter) modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--filter-name \"QD\" --filter \"QD < $params.qd_filter\" "])

// Initialize varaint annotation associated channels
def tools               = params.annotate_tools    ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
def snpeff_db           = params.snpeff_db         ?: Channel.empty()
def vep_cache_version   = params.vep_cache_version ?: Channel.empty()
def vep_genome          = params.vep_genome        ?: Channel.empty()
def vep_species         = params.vep_species       ?: Channel.empty()
def snpeff_cache        = params.snpeff_cache      ? params.snpeff_cache : []
def vep_cache           = params.vep_cache         ? params.vep_cache : []

// Include all the modules required for the pipeline
include { CAT_FASTQ }               from '../modules/nf-core/modules/cat/fastq/main'               addParams(options: modules['cat_fastq'])
include { GATK4_BASERECALIBRATOR }  from '../modules/nf-core/modules/gatk4/baserecalibrator/main'  addParams(options: modules['gatk_baserecalibrator'])
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/modules/gatk4/bedtointervallist/main' addParams(options: modules['gatk_bedtointervallist'])
include { GATK4_HAPLOTYPECALLER   } from '../modules/local/gatk4/haplotypecaller/main'             addParams(options: modules['gatk_haplotypecaller'])
include { GATK4_INTERVALLISTTOOLS } from '../modules/nf-core/modules/gatk4/intervallisttools/main' addParams(options: modules['gatk_intervallisttools'])
include { GATK4_MERGEVCFS         } from '../modules/nf-core/modules/gatk4/mergevcfs/main'         addParams(options: modules['gatk_mergevcfs'])
include { GATK4_INDEXFEATUREFILE  } from '../modules/nf-core/modules/gatk4/indexfeaturefile/main'  addParams(options: modules['gatk_indexfeaturefile'])
include { GATK4_VARIANTFILTRATION } from '../modules/nf-core/modules/gatk4/variantfiltration/main' addParams(options: modules['gatk_variantfilter'])
include { SAMTOOLS_INDEX }          from '../modules/nf-core/modules/samtools/index/main'          addParams(options: modules['samtools_index_genome'])

// Include all the subworkflows required for the pipeline

// Build the genome index and other reference files
include { PREPARE_GENOME }          from '../subworkflows/local/prepare_genome'                    addParams(
    genome_options:     publish_genome_options,
    index_options:      publish_index_options,
    star_index_options: modules['star_genomegenerate'],
    star_untar_options: untar_options
)
// Align reads to genome and sort and index the alignment file
include { ALIGN_STAR }              from '../subworkflows/nf-core/align_star'                      addParams(
    align_options: modules['star_align'],
    samtools_index_options: modules['samtools_index_genome'],
    samtools_sort_options:  modules['samtools_sort_genome'],
    samtools_stats_options: modules['samtools_index_genome'],
    seq_platform: params.seq_platform
)
// Mark duplicates in the BAM file
include { MARKDUPLICATES }          from '../subworkflows/nf-core/markduplicates'                  addParams(
    markduplicates_options: modules['picard_markduplicates'],
    samtools_index_options: modules['picard_markduplicates_samtools'],
    samtools_stats_options: modules['picard_markduplicates_samtools']
)
// Subworkflow - splits reads that contain Ns in their cigar string
include { SPLITNCIGAR } from '../subworkflows/local/splitncigar'                                    addParams(
    gatk_splitncigar_options: modules['gatk_splitncigar_options'],
    samtools_index_options: modules['samtools_index_genome'],
    samtools_merge_options: modules['samtools_merge_genome']
)

// Estimate and correct systematic bias
include { RECALIBRATE }             from '../subworkflows/nf-core/recalibrate'                     addParams(
    applybqsr_options: modules['gatk_applybqsr'],
    samtools_index_options: modules['samtools_index_recalibrate'],
    samtools_stats_options: modules['samtools_stats_recalibrate']
)
// Annotate variants using snpEff or VEP or both
include { ANNOTATE } from '../subworkflows/local/annotate' addParams(
    annotation_cache:                  params.annotation_cache,
    bgziptabix_merge_vep_options:      modules['bgziptabix_merge_vep'],
    bgziptabix_snpeff_options:         modules['bgziptabix_snpeff'],
    bgziptabix_vep_options:            modules['bgziptabix_vep'],
    merge_vep_options:                 modules['merge_vep'],
    snpeff_options:                    modules['snpeff'],
    snpeff_tag:                        "${modules['snpeff'].tag_base}.${params.genome}",
    vep_options:                       modules['ensemblvep'],
    vep_tag:                           "${modules['ensemblvep'].tag_base}.${params.genome}"
)

// Info required for completion email and summary
def multiqc_report = []

workflow RNAVAR {

    ch_versions = Channel.empty()
    ch_fastq    = Channel.empty()

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
            meta.id = meta.id.split('_')[0..-2].join('_')
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

    // PREPARE THE INTERVAL LIST FROM GTF FILE
    ch_interval_list = Channel.empty()
    GATK4_BEDTOINTERVALLIST(ch_genome_bed, PREPARE_GENOME.out.dict)
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first().ifEmpty(null))

    // MODULE: IntervalListTools from GATK4
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
            PREPARE_GENOME.out.gtf
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
        bam_splitncigar         = Channel.empty()
        SPLITNCIGAR(ch_genome_bam, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fai, PREPARE_GENOME.out.dict, ch_interval_list_split)
        bam_splitncigar         = SPLITNCIGAR.out.bam.join(SPLITNCIGAR.out.bai, by: [0])
        ch_versions             = ch_versions.mix(SPLITNCIGAR.out.versions.first().ifEmpty(null))

        // MODULE: BaseRecalibrator from GATK4
        ch_bqsr_table = Channel.empty()
        known_sites     = Channel.from([params.dbsnp, params.known_indels]).collect()
        known_sites_tbi = Channel.from([params.dbsnp_tbi, params.known_indels_tbi]).collect()

        ch_interval_list_baserecalibrator = ch_interval_list.map{ meta, bed -> [bed] }.collect()

        GATK4_BASERECALIBRATOR(
            bam_splitncigar,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict,
            ch_interval_list_baserecalibrator,
            known_sites,
            known_sites_tbi
        )
        ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table
        ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first().ifEmpty(null))

        // MODULE: ApplyBaseRecalibrator from GATK4
        bam_applybqsr       = bam_splitncigar.join(ch_bqsr_table, by: [0])
        bam_recalibrated    = Channel.empty()
        bam_recalibrated_qc = Channel.empty()

        ch_interval_list_applybqsr = ch_interval_list.map{ meta, bed -> [bed] }.collect()

        RECALIBRATE(
            ('bamqc' in params.skip_qc),
            ('samtools' in params.skip_qc),
            bam_applybqsr,
            PREPARE_GENOME.out.dict,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.fasta,
            ch_interval_list_applybqsr
        )

        bam_recalibrated    = RECALIBRATE.out.bam
        bam_recalibrated_qc = RECALIBRATE.out.qc
        ch_versions         = ch_versions.mix(RECALIBRATE.out.versions.first().ifEmpty(null))

        // MODULE: HaplotypeCaller from GATK4
        interval_flag = params.no_intervals
        haplotypecaller_vcf = Channel.empty()

        haplotypecaller_interval_bam = bam_recalibrated.combine(ch_interval_list_split)
            .map{ meta, bam, bai, interval_list ->
                new_meta = meta.clone()
                new_meta.id = meta.id
                [new_meta, bam, bai, interval_list]
            }

        GATK4_HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            params.dbsnp,
            params.dbsnp_tbi,
            PREPARE_GENOME.out.dict,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            interval_flag
        )

        haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
            .map{ meta, vcf ->
                meta.id = meta.id
                [meta, vcf]}
            .groupTuple()

        ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first().ifEmpty(null))
        use_ref_dict = true

        GATK4_MERGEVCFS(
            haplotypecaller_raw,
            PREPARE_GENOME.out.dict,
            use_ref_dict
        )
        haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
        ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first().ifEmpty(null))

        GATK4_INDEXFEATUREFILE(
            haplotypecaller_vcf
        )
        haplotypecaller_vcf_tbi = haplotypecaller_vcf.join(GATK4_INDEXFEATUREFILE.out.index, by: [0])
        ch_versions  = ch_versions.mix(GATK4_INDEXFEATUREFILE.out.versions.first().ifEmpty(null))
        final_vcf_tbi = haplotypecaller_vcf_tbi

        // MODULE: VariantFiltration from GATK4
        if (!params.skip_variantfiltration) {

            GATK4_VARIANTFILTRATION(
                haplotypecaller_vcf_tbi,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict
            )

            filtered_vcf_tbi = GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi, by: [0])
            final_vcf_tbi = filtered_vcf_tbi
            ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first().ifEmpty(null))
        }

        if(!params.skip_variantannotation && ('merge' in tools || 'snpeff' in tools || 'vep' in tools)) {
            ANNOTATE(
                final_vcf_tbi,
                tools,
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
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
