#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnavar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnavar
    Website: https://nf-co.re/rnavar
    Slack  : https://nfcore.slack.com/channels/rnavar
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.bcftools_annotations     = getGenomeAttribute('bcftools_annotations')
params.bcftools_annotations_tbi = getGenomeAttribute('bcftools_annotations_tbi')
params.dbsnp                    = getGenomeAttribute('dbsnp')
params.dbsnp_tbi                = getGenomeAttribute('dbsnp_tbi')
params.dict                     = getGenomeAttribute('dict')
params.exon_bed                 = getGenomeAttribute('exon_bed')
params.fasta                    = getGenomeAttribute('fasta')
params.fasta_fai                = getGenomeAttribute('fasta_fai')
params.gff                      = getGenomeAttribute('gff')
params.gtf                      = getGenomeAttribute('gtf')
params.known_indels             = getGenomeAttribute('known_indels')
params.known_indels_tbi         = getGenomeAttribute('known_indels_tbi')
params.snpeff_db                = getGenomeAttribute('snpeff_db')
params.star_index               = getGenomeAttribute('star')
params.vep_cache_version        = getGenomeAttribute('vep_cache_version')
params.vep_genome               = getGenomeAttribute('vep_genome')
params.vep_species              = getGenomeAttribute('vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAVAR                          } from './workflows/rnavar'
include { ANNOTATION_CACHE_INITIALISATION } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'
include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { PIPELINE_COMPLETION             } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { PREPARE_GENOME                  } from './subworkflows/local/prepare_genome'

// MULTIQC
include { MULTIQC                         } from './modules/nf-core/multiqc'
include { getWorkflowVersion              } from 'plugin/nf-core-utils'
include { processVersionsFromFile         } from 'plugin/nf-core-utils'
include { methodsDescriptionText          } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { paramsSummaryMap                } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RNAVAR {
    take:
    samplesheet
    align

    main:
    reports = Channel.empty()
    versions = Channel.empty()

    if (params.gtf && params.gff) {
        error("Using both --gtf and --gff is not supported. Please use only one of these parameters")
    }
    else if (!params.gtf && !params.gff) {
        error("Missing required parameters: --gtf or --gff")
    }

    if (params.extract_umi && !params.umitools_bc_pattern) {
        error("Expected --umitools_bc_pattern when --extract_umi is specified.")
    }

    if (!params.skip_baserecalibration && !params.dbsnp && !params.known_indels) {
        error("Known sites are required for performing base recalibration. Supply them with either --dbsnp and/or --known_indels or disable base recalibration with --skip_baserecalibration")
    }

    PREPARE_GENOME(
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.dbsnp,
        params.dbsnp_tbi,
        params.dict,
        params.exon_bed,
        params.fasta,
        params.fasta_fai,
        params.gff,
        params.gtf,
        params.known_indels,
        params.known_indels_tbi,
        params.star_index,
        params.feature_type,
        params.skip_exon_bed_check,
        align,
    )

    versions = versions.mix(PREPARE_GENOME.out.versions)

    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, rnavar will download the cache
        ensemblvep_info = Channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, params.vep_species, params.vep_cache_version])
        snpeff_info = Channel.of([[id: "${params.snpeff_db}"], params.snpeff_db])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map { _meta, cache -> [cache] }

        versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    }
    else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_db,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_custom_args,
            "Please refer to https://nf-co.re/rnavar/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
        )

        snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
        vep_cache = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    vep_extra_files = []

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
        vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }
    else if (params.dbnsfp && !params.dbnsfp_tbi) {
        System.err.println("DBNSFP: ${params.dbnsfp} has been provided with `--dbnsfp, but no dbnsfp_tbi has")
        System.err.println("cf: https://nf-co.re/rnavar/parameters/#dbnsfp")
        error("Execution halted due to dbnsfp inconsistency.")
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    //
    // WORKFLOW: Run pipeline
    //
    RNAVAR(
        samplesheet,
        PREPARE_GENOME.out.bcfann,
        PREPARE_GENOME.out.bcfann_tbi,
        params.bcftools_header_lines ? Channel.fromPath(params.bcftools_header_lines).collect() : Channel.empty(),
        PREPARE_GENOME.out.dbsnp,
        PREPARE_GENOME.out.dbsnp_tbi,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.exon_bed,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fasta_fai,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.known_sites,
        PREPARE_GENOME.out.known_sites_tbi,
        PREPARE_GENOME.out.star_index,
        snpeff_cache,
        params.snpeff_db,
        params.vep_genome,
        params.vep_species,
        params.vep_cache_version,
        params.vep_include_fasta,
        vep_cache,
        vep_extra_files,
        params.seq_center ?: [],
        params.seq_platform ?: [],
        params.aligner,
        params.bam_csi_index,
        params.extract_umi,
        params.generate_gvcf,
        params.skip_multiqc,
        params.skip_baserecalibration,
        params.skip_intervallisttools,
        params.skip_variantannotation,
        params.skip_variantfiltration,
        params.star_ignore_sjdbgtf,
        params.tools ?: "no_tools",
    )

    reports = reports.mix(RNAVAR.out.reports)
    versions = versions.mix(RNAVAR.out.versions)

    emit:
    reports  // channel: qc reports for multiQC
    versions // channel: [ path(versions.yml) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNAVAR(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.align,
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel
        .topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def collated_versions = softwareVersionsToYAML(NFCORE_RNAVAR.out.versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'rnavar_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    def multiqc_report = Channel.empty()

    // MULTIQC
    if (!params.skip_multiqc) {
        def multiqc_files = Channel.empty()

        def multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        def multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        def multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        def multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def methods_description = Channel.value(methodsDescriptionText(multiqc_custom_methods_description))

        multiqc_files = multiqc_files.mix(NFCORE_RNAVAR.out.reports)
        multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(collated_versions)
        multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC(
            multiqc_files.collect(),
            multiqc_config.toList(),
            multiqc_custom_config.toList(),
            multiqc_logo.toList(),
            [],
            [],
        )
        multiqc_report = MULTIQC.out.report.toList()
    }


    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions.unique().map { version -> processVersionsFromYAML(version) }.unique().mix(Channel.of(workflowVersionToYAML()))
}

//
// Process versions from topic channel
//
def topicVersionToYAML(taskProcess, tools, versions) {
    def toolsVersions = [tools, versions]
        .transpose()
        .collect { k, v -> "${k}: ${v}" }
    return """
    |${taskProcess.tokenize(':').last()}:
    |  ${toolsVersions.join('\n|  ')}
    """.stripMargin().trim()
}

//
// Get workflow summary for MultiQC
//
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params
        .keySet()
        .each { group ->
            def group_params = summary_params.get(group)
            // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>${group}</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                group_params
                    .keySet()
                    .sort()
                    .each { param ->
                        summary_section += "        <dt>${param}</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                    }
                summary_section += "    </dl>\n"
            }
        }

    def yaml_file_text = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n" as String
    yaml_file_text += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text += "plot_type: 'html'\n"
    yaml_file_text += "data: |\n"
    yaml_file_text += "${summary_section}"

    return yaml_file_text
}
