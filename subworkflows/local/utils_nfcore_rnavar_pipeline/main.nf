//
// Subworkflow with functionality specific to the nf-core/rnavar pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from 'plugin/nf-core-utils'
include { completionSummary       } from 'plugin/nf-core-utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir // string: The output directory where the results will be saved
    input // string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message
    bam_csi_index // params: params.bam_csi_index
    dbsnp // params: params.dbsnp
    gff // params: params.gff
    gtf // params: params.gtf
    known_indels // params: params.known_indels
    tools // array: list of tools to run, determined by params on launch time
    umitools_bc_pattern // params: params.umitools_bc_pattern

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                    \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                    \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/rnavar ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/rnavar/blob/master/CITATIONS.md
"""

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Extra informations
    //
    def extra_text = """
\033[1;37mExtra informations\033[0m
\033[0;34m  Tools selected to be run  :\033[0;32m ${tools.join(",")} \033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    if (params.monochrome_logs) {
        extra_text = extra_text.replaceAll(/\033\[[0-9;]*m/, '')
    }
    log.info(extra_text)

    //
    // Fails for missing params
    //
    if (gtf && gff) {
        error("Using both --gtf and --gff is not supported. Please use only one of these parameters")
    }
    else if (!gtf && !gff) {
        error("Missing required parameters: --gtf or --gff")
    }

    if (('umitools' in tools) && !umitools_bc_pattern) {
        error("Expected --umitools_bc_pattern when --tools umitools is specified.")
    }

    if (('baserecalibrator' in tools) && !dbsnp && !known_indels) {
        error("Known sites are required for performing base recalibration. Supply them with either --dbsnp and/or --known_indels or disable base recalibration with --skip_baserecalibration")
    }

    if ('variantfiltration' in tools && bam_csi_index) {
        log.warn('GATK4_VARIANTFILTRATION will not run with params.bam_csi_index')
    }

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    def samplesheetList = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    def bool_align = samplesheetList.find { _meta, fastq_1, _fastq_2, _bam, _bai, _cram, _crai, _vcf, _tbi ->
        fastq_1
    }
        ? true
        : false

    ch_samplesheet = channel.fromList(samplesheetList)
        .map { meta, fastq_1, fastq_2, bam, bai, cram, crai, vcf, tbi ->
            def new_meta = meta + [single_end: !fastq_2]
            [meta.id, new_meta, fastq_1, fastq_2, bam, bai, cram, crai, vcf, tbi]
        }

    emit:
    align       = bool_align
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email // string: email address
    email_on_fail // string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir // path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report // string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Function to check samples are internally consistent after being grouped
//
def checkSamplesAfterGrouping(input) {
    def (_ids, metas, fastqs_1, fastqs_2, bams, bais, crams, crais, vcfs, tbis) = input

    def fastqs_1_list = fastqs_1.findAll { it -> it != [] }
    def fastqs_2_list = fastqs_2.findAll { it -> it != [] }
    def bam_list = bams.findAll { it -> it != [] }
    def bai_list = bais.findAll { it -> it != [] }
    def cram_list = crams.findAll { it -> it != [] }
    def crai_list = crais.findAll { it -> it != [] }
    def vcf_list = vcfs.findAll { it -> it != [] }
    def tbi_list = tbis.findAll { it -> it != [] }

    def alignment_file_list = bam_list + cram_list
    if (alignment_file_list.size > 1) {
        error("Please check input samplesheet -> Multiple BAM/CRAM files per sample are not supported: ${metas[0].id}")
    }

    // Preserve paired-end structure by interleaving fastq_1 and fastq_2 files
    def fastqs = []
    if (fastqs_1_list.size() > 0 && fastqs_2_list.size() > 0) {
        // Paired-end: interleave the files to maintain pairing
        def max_size = Math.max(fastqs_1_list.size(), fastqs_2_list.size())
        (0..<max_size).each { i ->
            if (i < fastqs_1_list.size()) {
                fastqs.add(fastqs_1_list[i])
            }
            if (i < fastqs_2_list.size()) {
                fastqs.add(fastqs_2_list[i])
            }
        }
    }
    else if (fastqs_1_list.size() > 0) {
        // Single-end: just use fastq_1 files
        fastqs = fastqs_1_list
    }

    if (alignment_file_list.size() == 1 && fastqs.size() > 0) {
        error("Please check input samplesheet -> Detected FASTQ and BAM/CRAM files for the same sample. Please provide only one file type per sample: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs, bam_list[0] ?: [], bai_list[0] ?: [], cram_list[0] ?: [], crai_list[0] ?: [], vcf_list[0] ?: [], tbi_list[0] ?: []]
}

// Define list of tools to run
def defineToolsList(bam_csi_index, input_skip, input_tools) {

    // opt in tools
    def tools_list = input_tools ? input_tools.tokenize(',') : []

    // opt out tools
    def skip_list = input_skip ? input_skip.tokenize(',') : []

    // Specific tools not to execute depending of params
    // No variantfiltration if bam_csi_index as GATK4_VARIANTFILTRATION does not support csi index
    if (bam_csi_index) {
        tools_list = tools_list - 'variantfiltration'
    }

    return (tools_list - skip_list).sort().unique()
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
