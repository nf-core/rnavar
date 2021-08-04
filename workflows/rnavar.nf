/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnavar.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.gtf,
    params.gff,
    params.gene_bed,
    params.star_index]

for (param in checkPathParamList) {if (param) file(param, checkIfExists: true)}

// Check mandatory parameters
if (params.input) ch_input = file(params.input)
else exit 1, 'Input samplesheet not specified!'

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
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams(options: [publish_files : ['tsv':'']])

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams(options: [:])

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

modules['multiqc'].args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams(options: modules['fastqc'])
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams(options: modules['multiqc'])

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices = params.aligner }

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

if (!params.save_reference) modules['star_genomegenerate']['publish_files'] = false

modules['samtools_index_genome'].args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

if (['star'].contains(params.aligner)) {
    if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
        modules['samtools_sort_genome'].publish_files.put('bam','')
        modules['samtools_index_genome'].publish_files.put('bai','')
        modules['samtools_index_genome'].publish_files.put('csi','')
    }
} else {
    if (params.save_align_intermeds || params.skip_markduplicates) {
        modules['samtools_sort_genome'].publish_files.put('bam','')
        modules['samtools_index_genome'].publish_files.put('bai','')
        modules['samtools_index_genome'].publish_files.put('csi','')
    }
}

def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type

if (!params.save_merged_fastq) modules['cat_fastq']['publish_files'] = false

modules['star_align'].args += params.save_unaligned ? Utils.joinModuleArgs(['--outReadsUnmapped Fastx']) : ''
if (params.save_align_intermeds) modules['star_align'].publish_files.put('bam','')
if (params.save_unaligned)       modules['star_align'].publish_files.put('fastq.gz','unmapped')

modules['picard_markduplicates_samtools'].args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''
modules['gatk_intervallisttools'].args += Utils.joinModuleArgs(["--SCATTER_COUNT $params.scatter_count"])

if (params.window)    modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--window $params.window"])
if (params.cluster)   modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--cluster $params.cluster"])
if (params.fs_filter) modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--filter-name \"FS\" --filter \"FS > $params.fs_filter\" "])
if (params.qd_filter) modules['gatk_variantfilter'].args += Utils.joinModuleArgs(["--filter-name \"QD\" --filter \"QD < $params.qd_filter\" "])

include { CAT_FASTQ }               from '../modules/nf-core/modules/cat/fastq/main'               addParams(options: modules['cat_fastq'])
include { GATK4_BASERECALIBRATOR }  from '../modules/nf-core/modules/gatk4/baserecalibrator/main'  addParams(options: modules['gatk_baserecalibrator'])
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/modules/gatk4/bedtointervallist/main' addParams(options: [:])
include { GATK4_HAPLOTYPECALLER   } from '../modules/local/gatk4/haplotypecaller/main'             addParams(options: modules['gatk_haplotypecaller'])
include { GATK4_INTERVALLISTTOOLS } from '../modules/nf-core/modules/gatk4/intervallisttools/main' addParams(options: modules['gatk_intervallisttools'])
include { GATK4_MERGEVCFS         } from '../modules/local/gatk4/mergevcfs/main'                   addParams(options: modules['gatk_mergevcfs'])
include { GATK4_SPLITNCIGARREADS }  from '../modules/local/gatk4/splitncigarreads/main'            addParams(options: modules['gatk_splitncigar_options'])
include { GATK4_VARIANTFILTRATION } from '../modules/local/gatk4/variantfiltration/main'           addParams(options: modules['gatk_variantfilter'])
include { SAMTOOLS_INDEX }          from '../modules/nf-core/modules/samtools/index/main'          addParams(options: modules['samtools_index_genome'])

include { PREPARE_GENOME }          from '../subworkflows/local/prepare_genome'                    addParams(
    genome_options:     publish_genome_options,
    index_options:      publish_index_options,
    star_index_options: modules['star_genomegenerate']
)
include { ALIGN_STAR }              from '../subworkflows/nf-core/align_star'                      addParams(
    align_options: modules['star_align'],
    samtools_index_options: modules['samtools_index_genome'],
    samtools_sort_options:  modules['samtools_sort_genome'],
    samtools_stats_options: modules['samtools_index_genome'],
    seq_platform: params.seq_platform
)
include { MARKDUPLICATES }          from '../subworkflows/nf-core/markduplicates'                  addParams(
    markduplicates_options: modules['picard_markduplicates'],
    samtools_index_options: modules['picard_markduplicates_samtools'],
    samtools_stats_options: modules['picard_markduplicates_samtools']
)
include { RECALIBRATE }             from '../subworkflows/nf-core/recalibrate'                     addParams(
    applybqsr_options:      modules['applybqsr'],
    merge_bam_options:      modules['merge_bam_recalibrate'],
    qualimap_bamqc_options: modules['qualimap_bamqc_recalibrate'],
    samtools_index_options: modules['samtools_index_recalibrate'],
    samtools_stats_options: modules['samtools_stats_recalibrate']
)

// Info required for completion email and summary
def multiqc_report = []

workflow RNAVAR {

    ch_software_versions = Channel.empty()

    PREPARE_GENOME(prepareToolIndices, biotype)
    ch_genome_gtf = Channel.from([id:'genome.bed']).combine(PREPARE_GENOME.out.gene_bed)

    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    INPUT_CHECK(ch_input)
        .map{ meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [meta, fastq]}
        .groupTuple(by: [0])
        .branch{ meta, fastq ->
            single  : fastq.size() == 1
                return [meta, fastq.flatten()]
            multiple: fastq.size() > 1
                return [meta, fastq.flatten()]}
        .set{ch_fastq}

    //
    // MODULE: Run FastQC
    //
    // FASTQC(ch_fastq)
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(ch_fastq.multiple)
        .mix(ch_fastq.single)
        .set{ch_cat_fastq}

    // PREPARE THE INTERVAL LIST
    ch_interval_list = Channel.empty()
    GATK4_BEDTOINTERVALLIST(ch_genome_gtf, PREPARE_GENOME.out.dict)
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_software_versions = ch_software_versions.mix(GATK4_BEDTOINTERVALLIST.out.version.first().ifEmpty(null))

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

    if (!params.skip_alignment && params.aligner == 'star') {
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
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))

        // SUBWORKFLOW: Mark duplicates with Picard
        ch_markduplicates_multiqc = Channel.empty()
        if (!params.skip_alignment && !params.skip_markduplicates) {
            MARKDUPLICATES(ch_genome_bam)
            ch_genome_bam             = MARKDUPLICATES.out.bam
            ch_genome_bam_index       = MARKDUPLICATES.out.bai
            ch_samtools_stats         = MARKDUPLICATES.out.stats
            ch_samtools_flagstat      = MARKDUPLICATES.out.flagstat
            ch_samtools_idxstats      = MARKDUPLICATES.out.idxstats
            ch_markduplicates_multiqc = MARKDUPLICATES.out.metrics
            if (params.bam_csi_index) ch_genome_bam_index = MARKDUPLICATES.out.csi
            ch_software_versions      = ch_software_versions.mix(MARKDUPLICATES.out.picard_version.first().ifEmpty(null))
        }
        // MODULE: SplitNCigarReads from GATK4
        // Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
        bam_splitncigar = Channel.empty()
        if (!params.skip_splitncigar) {
            GATK4_SPLITNCIGARREADS(ch_genome_bam, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fai, PREPARE_GENOME.out.dict)
            bam_splitncigar      = GATK4_SPLITNCIGARREADS.out.bam
            ch_software_versions = ch_software_versions.mix(GATK4_SPLITNCIGARREADS.out.version.first().ifEmpty(null))
        }

        ch_bqsr_table = Channel.empty()

        if (!params.skip_basecalibrator) {
            known_sites     = Channel.from([params.dbsnp_vcf, params.known_indels]).collect()
            known_sites_tbi = Channel.from([params.dbsnp_vcf_index, params.known_indels_index]).collect()

            ch_interval_list_baserecalibrator = ch_interval_list.map{ meta, bed -> [bed] }

            GATK4_BASERECALIBRATOR(
                bam_splitncigar,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict,
                ch_interval_list_baserecalibrator,
                known_sites,
                known_sites_tbi
            )
            ch_bqsr_table = GATK4_BASERECALIBRATOR.out.table
        }

        bam_applybqsr = bam_splitncigar.join(ch_bqsr_table, by: [0])

        bam_recalibrated    = Channel.empty()
        bam_recalibrated_qc = Channel.empty()

        if (!params.skip_applybqsr) {

            ch_interval_list_applybqsr = ch_interval_list.map{ meta, bed -> [bed] }
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

        }

        ch_interval_list_split = Channel.empty()

        if (!params.skip_intervallisttools) {
            GATK4_INTERVALLISTTOOLS(ch_interval_list)
            ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ meta, bed -> [bed] }.flatten()
        }
        else ch_interval_list_split = ch_interval_list

        interval_flag = false

        haplotypecaller_vcf = Channel.empty()

        if (!params.skip_haplotypecaller) {

            haplotypecaller_interval_bam = bam_recalibrated.combine(ch_interval_list_split)
                .map{ meta, bam, bai, interval_list ->
                    new_meta = meta.clone()
                    new_meta.id = meta.id + "_" + interval_list.baseName
                    [new_meta, bam, bai, interval_list]}

            GATK4_HAPLOTYPECALLER(
                haplotypecaller_interval_bam,
                params.dbsnp_vcf,
                params.dbsnp_vcf_index,
                PREPARE_GENOME.out.dict,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                interval_flag
            )

            haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
                .map{ meta,vcf ->
                    meta.id = meta.sample
                    [meta, vcf]}
                .groupTuple()

            use_ref_dict = true
            GATK4_MERGEVCFS(
                haplotypecaller_raw,
                PREPARE_GENOME.out.dict,
                use_ref_dict
            )

            haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
        }

        if (!params.skip_variantfiltration) {

            GATK4_VARIANTFILTRATION(
                haplotypecaller_vcf,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict
            )

            filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
        }

    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map{it -> if (it) [it.baseName, it]}
        .groupTuple()
        .map{it[1][0]}
        .flatten()
        .collect()
        .set{ch_software_versions}

    GET_SOFTWARE_VERSIONS(ch_software_versions.map{it}.collect())

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRnavar.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (ch_multiqc_files.collect())
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
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
