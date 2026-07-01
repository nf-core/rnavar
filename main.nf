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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAVAR                           } from './workflows/rnavar'
include { PREPARE_GENOME                   } from './subworkflows/local/prepare_genome'
include { CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF } from './subworkflows/nf-core/cache_download_ensemblvep_snpeff'
include { UTILS_ANNOTATION_CACHE           } from './subworkflows/nf-core/utils_annotation_cache'
include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_rnavar_pipeline'

// MULTIQC
include { MULTIQC                          } from './modules/nf-core/multiqc'
include { methodsDescriptionText           } from './subworkflows/local/utils_nfcore_rnavar_pipeline'
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { softwareVersionsToYAML           } from 'plugin/nf-core-utils'

// tools selections
include { defineToolsList                  } from './subworkflows/local/utils_nfcore_rnavar_pipeline'

// references
include { getGenomeAttribute               } from 'plugin/nf-core-utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // Input options
    input: Path
    outdir: String
    tools: String?
    skip_tools: String?

    // References
    snpeff_cache: String = 's3://annotation-cache/snpeff_cache/'
    vep_cache: String = 's3://annotation-cache/vep_cache/'
    save_reference: Boolean = false
    feature_type: String = 'exon'

    // Sequence read information
    save_merged_fastq: Boolean = false
    read_length: Integer = 150

    // Alignment
    aligner: String = 'star'
    star_twopass: Boolean = true
    star_ignore_sjdbgtf: Boolean = false
    star_max_memory_bamsort: Integer = 0
    star_bins_bamsort: Integer = 50
    star_max_collapsed_junc: Integer = 1000000
    star_max_intron_size: Integer?
    seq_center: String?
    seq_platform: String = 'illumina'
    bam_csi_index: Boolean = false
    save_unaligned: Boolean = false
    save_align_intermeds: Boolean = false

    // Preprocessing of alignment
    remove_duplicates: Boolean = false
    umitools_extract_method: String = 'string'
    umitools_bc_pattern: String?
    umitools_bc_pattern2: String?
    umitools_umi_separator: String?

    // Variant calling
    no_intervals: Boolean = false

    // Variant annotation
    bcftools_annotations: String?
    bcftools_annotations_tbi: String?
    bcftools_columns: String?
    bcftools_header_lines: String?
    download_cache: Boolean = false
    dbnsfp: String?
    dbnsfp_consequence: String?
    dbnsfp_fields: String = "rs_dbSNP,HGVSc_VEP,HGVSp_VEP,1000Gp3_EAS_AF,1000Gp3_AMR_AF,LRT_score,GERP++_RS,gnomAD_exomes_AF"
    dbnsfp_tbi: String?
    outdir_cache: String?
    spliceai_indel: String?
    spliceai_indel_tbi: String?
    spliceai_snv: String?
    spliceai_snv_tbi: String?
    vep_custom_args: String = "--everything --filter_common --per_gene --total_length --offline --format vcf"
    vep_dbnsfp: Boolean?
    vep_include_fasta: Boolean = false
    vep_loftee: Boolean?
    vep_version: String = "115.2-1"
    vep_out_format: String = "vcf"
    vep_cache_preflight_check: Boolean = false
    vep_spliceai: Boolean?
    vep_spliceregion: Boolean?

    // GATK intervallist parameters
    gatk_interval_scatter_count: Integer = 25

    // GATK haplotypecaller parameters
    gatk_hc_call_conf: Integer = 20
    generate_gvcf: Boolean = false

    //GATK variant filter parameters
    gatk_vf_window_size: Integer = 35
    gatk_vf_cluster_size: Integer = 3
    gatk_vf_fs_filter: Float = 30.0
    gatk_vf_qd_filter: Float = 2.0

    // MultiQC options
    multiqc_config: Path?
    multiqc_title: String?
    multiqc_logo: Path?
    multiqc_methods_description: String?
    max_multiqc_email_size: String = '25.MB'

    // Boilerplate options
    email: String?
    email_on_fail: String?
    plaintext_email: Boolean = false
    help: Boolean = false
    help_full: Boolean = false
    show_hidden: Boolean = false
    version: Boolean = false
    modules_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'

    // Config options
    config_profile_name: String?
    config_profile_description: String?
    config_profile_contact: String?
    config_profile_url: String?

    // Schema validation default options
    validate_params: Boolean = true

    // GENOME PARAMETER VALUES
    genomes: Map = [:]
    dbsnp: Path? = getGenomeAttribute('dbsnp')
    dbsnp_tbi: Path? = getGenomeAttribute('dbsnp_tbi')
    dict: Path? = getGenomeAttribute('dict')
    exon_bed: Path? = getGenomeAttribute('exon_bed')
    fasta: Path = getGenomeAttribute('fasta')
    fasta_fai: Path? = getGenomeAttribute('fasta_fai')
    gff: Path? = getGenomeAttribute('gff')
    gtf: Path? = getGenomeAttribute('gtf')
    known_indels: Path? = getGenomeAttribute('known_indels')
    known_indels_tbi: Path? = getGenomeAttribute('known_indels_tbi')
    snpeff_db: String? = getGenomeAttribute('snpeff_db')
    star_index: Path? = getGenomeAttribute('star')
    vep_cache_version: Integer? = getGenomeAttribute('vep_cache_version')
    vep_genome: String? = getGenomeAttribute('vep_genome')
    vep_species: String? = getGenomeAttribute('vep_species')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    def tools = defineToolsList(
        params.bam_csi_index,
        params.skip_tools,
        params.tools,
    )

    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.bam_csi_index,
        params.dbsnp,
        params.gff,
        params.gtf,
        params.known_indels,
        tools,
        params.umitools_bc_pattern,
    )


    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, rnavar will download the cache
        ensemblvep_info = channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, params.vep_species, params.vep_cache_version])
        snpeff_info = channel.of([[id: "${params.snpeff_db}"], params.snpeff_db])
        CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF(ensemblvep_info, snpeff_info, params.vep_cache_preflight_check)
        snpeff_cache = CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF.out.snpeff_cache
        vep_cache = CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF.out.ensemblvep_cache
    }
    else {
        // Looks for cache information either locally or on the cloud
        UTILS_ANNOTATION_CACHE(
            params.vep_cache,
            params.vep_cache_version,
            params.vep_custom_args,
            params.vep_genome,
            params.vep_species,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_db,
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            "Please refer to https://nf-co.re/rnavar/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
        )

        snpeff_cache = UTILS_ANNOTATION_CACHE.out.snpeff_cache
        vep_cache = UTILS_ANNOTATION_CACHE.out.ensemblvep_cache
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

    // WORKFLOW: Run main workflow
    NFCORE_RNAVAR(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.align,
        snpeff_cache,
        vep_cache,
        vep_extra_files,
        tools,
    )

    def collated_versions = softwareVersionsToYAML(
        softwareVersions: channel.topic("versions"),
        nextflowVersion: workflow.nextflow.version,
    ).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_' + 'rnavar_software_' + 'mqc_' + 'versions.yml',
        sort: true,
        newLine: true,
    )

    def collated_reports = channel.topic("multiqc_files")
        .map { _meta, _process, _tool, reports -> reports }

    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    def multiqc_report = channel.empty()

    // MULTIQC
    def multiqc_files = channel.empty()

    multiqc_files = multiqc_files.mix(collated_versions)
    multiqc_files = multiqc_files.mix(collated_reports)

    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    def multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def methods_description = channel.value(methodsDescriptionText(multiqc_custom_methods_description))

    multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        multiqc_files.flatten().collect().map { files ->
            [
                [id: 'rnavar'],
                files,
                params.multiqc_config
                    ? file(params.multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }.filter { ('multiqc' in tools) }
    )
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        multiqc_report,
    )

    publish:
    multiqc = MULTIQC.out.data.mix(MULTIQC.out.plots, MULTIQC.out.report)
    reports = channel.topic("multiqc_files").filter { _meta, _process, tool, _file ->
        return !(tool == 'gatk4' || (tool == 'snpeff' && !('snpeff' in tools)))
    }
}

output {
    multiqc {
        path "reports/multiqc"
    }
    reports {
        path { meta, _process, tool, file ->
            file >> "reports/${tool}/${meta.id}/"
        }
    }
}

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
    snpeff_cache
    vep_cache
    vep_extra_files
    tools

    main:
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
        align,
        params.genome ?: "genome",
        tools,
    )

    // WORKFLOW: Run pipeline
    RNAVAR(
        samplesheet,
        PREPARE_GENOME.out.bcfann,
        PREPARE_GENOME.out.bcfann_tbi,
        params.bcftools_columns ? channel.fromPath(params.bcftools_columns).collect() : false,
        params.bcftools_header_lines ? channel.fromPath(params.bcftools_header_lines).collect() : channel.empty(),
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
        params.aligner,
        params.star_ignore_sjdbgtf,
        tools,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get workflow summary for MultiQC
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
