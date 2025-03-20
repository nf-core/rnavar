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

    main:

    ch_versions = Channel.empty()

    if(params.gtf && params.gff) {
        error("Using both --gtf and --gff is not supported. Please use only one of these parameters")
    } else if(!params.gtf && !params.gff) {
        error("Missing required parameters: --gtf or --gff")
    }

    if(params.extract_umi && !params.umitools_bc_pattern) {
        error("Expected --umitools_bc_pattern when --extract_umi is specified.")
    }

    if(!params.skip_baserecalibration && !params.dbsnp && !params.known_indels) {
        error("Known sites are required for performing base recalibration. Supply them with either --dbsnp and/or --known_sites or disable base recalibration with --skip_baserecalibration")
    }

    // Initialize fasta file with meta map:
    ch_fasta_raw      = params.fasta                   ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect()        : Channel.empty()

    // Initialize file channels based on params, defined in the params.genomes[params.genome] scope
    ch_dict_raw       = params.dict                    ? Channel.fromPath(params.dict).map{ it -> [ [id:it.baseName], it ] }.collect()         : Channel.empty()
    ch_fai_raw        = params.fasta_fai               ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:it.baseName], it ] }.collect()    : Channel.empty()
    ch_dbsnp_raw      = params.dbsnp                   ? Channel.fromPath(params.dbsnp).map { dbsnp -> [[id:dbsnp.baseName], dbsnp]}.collect() : Channel.value([])
    ch_known_indels_raw     = params.known_indels      ? Channel.fromPath(params.known_indels)                                                 : Channel.empty()
    ch_known_indels_tbi_raw = params.known_indels_tbi  ? Channel.fromPath(params.known_indels_tbi)                                             : Channel.empty()
    ch_gff            = params.gff                     ? Channel.fromPath(params.gff).map{ gff -> [ [ id:gff.baseName ], gff ] }.collect()     : Channel.empty()
    ch_gtf_raw        = params.gtf                     ? Channel.fromPath(params.gtf).map{ gtf -> [ [ id:gtf.baseName ], gtf ] }.collect()     : Channel.empty()
    ch_star_index_raw = params.star_index              ? Channel.fromPath(params.star_index).map { index -> [[id:index.baseName], index]}      : Channel.value([[],[]])
    ch_exon_bed_raw   = params.exon_bed                ? Channel.fromPath(params.exon_bed).map { it -> [[id:it.baseName], it]}                 : Channel.empty()

    // Initialize variant annotation associated channels
    snpeff_db         = params.snpeff_db          ?:  Channel.empty()
    vep_cache_version = params.vep_cache_version  ?:  Channel.empty()
    vep_genome        = params.vep_genome         ?:  Channel.empty()
    vep_species       = params.vep_species        ?:  Channel.empty()

    seq_platform      = params.seq_platform ?: []
    seq_center        = params.seq_center   ?: []

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
        ch_fasta_raw,
        ch_dict_raw,
        ch_fai_raw,
        ch_star_index_raw,
        ch_gff,
        ch_gtf_raw,
        ch_exon_bed_raw,
        ch_dbsnp_raw,
        ch_known_indels_raw,
        ch_known_indels_tbi_raw,
        params.feature_type)

    ch_fasta            = PREPARE_GENOME.out.fasta
    ch_star_index       = PREPARE_GENOME.out.star_index
    ch_gtf              = PREPARE_GENOME.out.gtf
    ch_dict             = PREPARE_GENOME.out.dict
    ch_fasta_fai        = PREPARE_GENOME.out.fasta_fai
    ch_exon_bed         = PREPARE_GENOME.out.exon_bed
    ch_dbsnp            = params.dbsnp && params.dbsnp.endsWith(".gz") ? ch_dbsnp_raw : PREPARE_GENOME.out.dbsnp
    ch_dbsnp_tbi        = params.dbsnp.toString().endsWith(".gz") && params.dbsnp_tbi
                                                ? Channel.fromPath(params.dbsnp_tbi).map { dbsnp -> [[id:dbsnp.baseName], dbsnp]}.collect()
                                                : PREPARE_GENOME.out.dbsnp_tbi
    ch_known_indels     = params.known_indels   ? PREPARE_GENOME.out.known_indels     : Channel.value([])
    ch_known_indels_tbi = params.known_indels   ? PREPARE_GENOME.out.known_indels_tbi : Channel.value([])

    vep_fasta = (params.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]

    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, sarek will download the cache
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        snpeff_info     = Channel.of([ [ id:"${params.snpeff_db}" ], params.snpeff_db ])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        ch_versions = ch_versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    } else {
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
            "Please refer to https://nf-co.re/rnavar/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

            snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
            vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    //
    // WORKFLOW: Run pipeline
    //
    RNAVAR(samplesheet,
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
        vep_fasta,
        vep_cache,
        vep_extra_files,
        seq_center,
        seq_platform)

    emit:
    multiqc_report = RNAVAR.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNAVAR(PIPELINE_INITIALISATION.out.samplesheet)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RNAVAR.out.multiqc_report
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
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
