//
// Subworkflow with functionality specific to the nf-core/rnavar pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { checkCondaChannels   } from 'plugin/nf-core-utils'
include { checkConfigProvided  } from 'plugin/nf-core-utils'
include { checkProfileProvided } from 'plugin/nf-core-utils'
include { completionEmail      } from 'plugin/nf-core-utils'
include { completionSummary    } from 'plugin/nf-core-utils'
include { dumpParametersToJSON } from 'plugin/nf-core-utils'
include { getWorkflowVersion   } from 'plugin/nf-core-utils'
include { imNotification       } from 'plugin/nf-core-utils'

include { paramsHelp           } from 'plugin/nf-schema'
include { paramsSummaryLog     } from 'plugin/nf-schema'
include { paramsSummaryMap     } from 'plugin/nf-schema'
include { samplesheetToList    } from 'plugin/nf-schema'
include { validateParameters   } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    // Print workflow version and exit on --version
    if (version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }

    // Dump pipeline parameters to a JSON file
    if (outdir) {
        dumpParametersToJSON(outdir, params)
    }

    // When running with Conda, warn if channels have not been set-up appropriately
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        checkCondaChannels()
    }

    checkConfigProvided()

    // Validate parameters and generate parameter summary to stdout
    //
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

    if (help || help_full) {
        help_options = [
            beforeText: before_text,
            afterText: after_text,
            command: command,
            showHidden: show_hidden,
            fullHelp: help_full,
        ]
        if (null) {
            help_options << [parametersSchema: null]
        }
        log.info(
            paramsHelp(
                help_options,
                params.help instanceof String && params.help != "true" ? params.help : "",
            )
        )
        exit(0)
    }

    checkProfileProvided(nextflow_cli_args)

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //

    summary_options = [:]
    if (null) {
        summary_options << [parametersSchema: null]
    }
    log.info(before_text)
    log.info(paramsSummaryLog(summary_options, workflow))
    log.info(after_text)

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if (validate_params) {
        validateOptions = [:]
        if (null) {
            validateOptions << [parametersSchema: null]
        }
        validateParameters(validateOptions)
    }

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // Create channel from input file provided through input
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
    samplesheet = ch_samplesheet
    align       = bool_align
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url //  string: hook URL for notifications
    multiqc_report //  string: Path to MultiQC report

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
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
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
    if (alignment_file_list.size() > 1) {
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
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
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

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
