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
params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.dict              = getGenomeAttribute('dict')
params.gtf               = getGenomeAttribute('gtf')
params.gff               = getGenomeAttribute('gff')
params.exon_bed          = getGenomeAttribute('exon_bed')
params.star_index        = getGenomeAttribute('star')
params.dbsnp             = getGenomeAttribute('dbsnp')
params.dbsnp_tbi         = getGenomeAttribute('dbsnp_tbi')
params.known_indels      = getGenomeAttribute('known_indels')
params.known_indels_tbi  = getGenomeAttribute('known_indels_tbi')
params.snpeff_db         = getGenomeAttribute('snpeff_db')
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')

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
include { methodsDescriptionText          } from './subworkflows/local/utils_nfcore_rnavar_pipeline'

// nf-core
include { paramsSummaryMultiqc            } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC                         } from './modules/nf-core/multiqc'

// plugin
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
        error("Known sites are required for performing base recalibration. Supply them with either --dbsnp and/or --known_sites or disable base recalibration with --skip_baserecalibration")
    }

    // Initialize file channels based on params
    ch_bcftools_annotations_raw = params.bcftools_annotations ? Channel.fromPath(params.bcftools_annotations) : Channel.empty()
    ch_bcftools_annotations_tbi_raw = params.bcftools_annotations_tbi ? Channel.fromPath(params.bcftools_annotations_tbi) : Channel.empty()
    ch_bcftools_header_lines = params.bcftools_header_lines ? Channel.fromPath(params.bcftools_header_lines).collect() : Channel.empty()
    ch_dbsnp_raw = params.dbsnp ? Channel.fromPath(params.dbsnp) : Channel.empty()
    ch_dbsnp_tbi_raw = params.dbsnp_tbi ? Channel.fromPath(params.dbsnp_tbi) : Channel.empty()
    ch_exon_bed_raw = params.exon_bed ? Channel.fromPath(params.exon_bed).map { it -> [[id: it.baseName], it] } : Channel.empty()
    ch_gff = params.gff ? Channel.fromPath(params.gff).map { gff -> [[id: gff.baseName], gff] }.collect() : Channel.empty()
    ch_gtf_raw = params.gtf ? Channel.fromPath(params.gtf).map { gtf -> [[id: gtf.baseName], gtf] }.collect() : Channel.empty()
    ch_known_indels_raw = params.known_indels ? Channel.fromPath(params.known_indels) : Channel.empty()
    ch_known_indels_tbi_raw = params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi) : Channel.empty()
    ch_star_index_raw = params.star_index ? Channel.fromPath(params.star_index).map { index -> [[id: index.baseName], index] } : Channel.value([[], []])

    seq_platform = params.seq_platform ?: []
    seq_center = params.seq_center ?: []

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

    PREPARE_GENOME(
        params.fasta,
        params.dict,
        params.fasta_fai,
        ch_star_index_raw,
        ch_gff,
        ch_gtf_raw,
        ch_exon_bed_raw,
        ch_bcftools_annotations_raw,
        ch_bcftools_annotations_tbi_raw,
        ch_dbsnp_raw,
        ch_dbsnp_tbi_raw,
        ch_known_indels_raw,
        ch_known_indels_tbi_raw,
        params.feature_type,
        params.skip_exon_bed_check,
        align,
    )

    ch_fasta = PREPARE_GENOME.out.fasta
    ch_star_index = PREPARE_GENOME.out.star_index
    ch_gtf = PREPARE_GENOME.out.gtf
    ch_dict = PREPARE_GENOME.out.dict
    ch_fasta_fai = PREPARE_GENOME.out.fasta_fai
    ch_exon_bed = PREPARE_GENOME.out.exon_bed
    ch_bcfann = params.bcftools_annotations ? PREPARE_GENOME.out.bcfann : Channel.value([])
    ch_bcfann_tbi = params.bcftools_annotations ? PREPARE_GENOME.out.bcfann_tbi : Channel.value([])
    ch_dbsnp = params.dbsnp ? PREPARE_GENOME.out.dbsnp : Channel.value([])
    ch_dbsnp_tbi = params.dbsnp ? PREPARE_GENOME.out.dbsnp_tbi : Channel.value([])
    ch_known_indels = params.known_indels ? PREPARE_GENOME.out.known_indels : Channel.value([])
    ch_known_indels_tbi = params.known_indels ? PREPARE_GENOME.out.known_indels_tbi : Channel.value([])

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

    //
    // WORKFLOW: Run pipeline
    //
    RNAVAR(
        samplesheet,
        ch_bcfann,
        ch_bcfann_tbi,
        ch_bcftools_header_lines,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_dict,
        ch_exon_bed,
        ch_fasta,
        ch_fasta_fai,
        ch_gtf,
        ch_known_indels,
        ch_known_indels_tbi,
        ch_star_index,
        snpeff_cache,
        params.snpeff_db,
        params.vep_genome,
        params.vep_species,
        params.vep_cache_version,
        params.vep_include_fasta,
        vep_cache,
        vep_extra_files,
        seq_center,
        seq_platform,
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
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNAVAR(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.align,
    )

    // Collate and save software versions
    def collated_versions = softwareVersionsToYAML(NFCORE_RNAVAR.out.versions).collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnavar_software_mqc_versions.yml', sort: true, newLine: true)

    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    def val_multiqc_report = Channel.empty()

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
        val_multiqc_report = MULTIQC.out.report.toList()
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
        val_multiqc_report,
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
