/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnavar.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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

if (!params.save_reference) modules['star_genomegenerate']['publish_files'] = false

modules['samtools_index_genome'].args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

if (params.save_align_intermeds) {
    modules['samtools_sort_genome'].publish_files.put('bam','')
    modules['samtools_index_genome'].publish_files.put('bai','')
    modules['samtools_index_genome'].publish_files.put('csi','')
}

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
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/modules/gatk4/bedtointervallist/main' addParams(options: modules['gatk_bedtointervallist'])
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
    samtools_index_options: modules['samtools_index_recalibrate'],
    samtools_stats_options: modules['samtools_stats_recalibrate']
)


// Info required for completion email and summary
def multiqc_report = []

workflow RNAVAR {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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
