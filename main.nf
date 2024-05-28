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

nextflow.enable.dsl = 2

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
params.snpeff_genome     = getGenomeAttribute('snpeff_genome')
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

    // Initialize fasta file with meta map:
    ch_fasta_raw      = params.fasta                   ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect()    : Channel.empty()

    // Initialize file channels based on params, defined in the params.genomes[params.genome] scope
    ch_dbsnp          = params.dbsnp                   ? Channel.fromPath(params.dbsnp).collect()                                          : Channel.value([])
    ch_known_indels   = params.known_indels            ? Channel.fromPath(params.known_indels).collect()                                   : Channel.value([])
    ch_gff            = params.gff                     ? Channel.fromPath(params.gff).collect()                                            : Channel.empty()
    ch_gtf_raw        = params.gtf                     ? Channel.fromPath(params.gtf).map{ gtf -> [ [ id:gtf.baseName ], gtf ] }.collect() : Channel.empty()

    // Initialize variant annotation associated channels
    snpeff_db         = params.snpeff_db          ?: Channel.empty()
    vep_cache_version = params.vep_cache_version  ?: Channel.empty()
    vep_genome        = params.vep_genome         ?: Channel.empty()
    vep_species       = params.vep_species        ?: Channel.empty()

    seq_platform      = params.seq_platform ?: []
    seq_center        = params.seq_center   ?: []

    // Initialize value channels based on params, defined in the params.genomes[params.genome] scope
    ascat_genome                = params.ascat_genome       ?:  Channel.empty()
    dbsnp_vqsr                  = params.dbsnp_vqsr         ?   Channel.value(params.dbsnp_vqsr)        : Channel.empty()
    known_indels_vqsr           = params.known_indels_vqsr  ?   Channel.value(params.known_indels_vqsr) : Channel.empty()
    known_snps_vqsr             = params.known_snps_vqsr    ?   Channel.value(params.known_snps_vqsr)   : Channel.empty()
    ngscheckmate_bed            = params.ngscheckmate_bed   ?   Channel.value(params.ngscheckmate_bed)  : Channel.empty()
    snpeff_db                   = params.snpeff_db          ?:  Channel.empty()
    vep_cache_version           = params.vep_cache_version  ?:  Channel.empty()
    vep_genome                  = params.vep_genome         ?:  Channel.empty()
    vep_species                 = params.vep_species        ?:  Channel.empty()

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
        ch_gff,
        ch_gtf_raw,
        feature_type)

    ch_dict             = params.dict           ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                                : PREPARE_GENOME.out.dict
    ch_fasta_fai        = params.fasta_fai      ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect()
                                                : PREPARE_GENOME.out.fasta_fai
    ch_exon_bed         = params.exon_bed       ? Channel.fromPath(params.exon_bed).map{ it -> [ [id:'exon_bed'], it ] }.collect()
                                                : PREPARE_GENOME.out.exon_bed
    ch_gtf              = params.gtf            ? Channel.fromPath(params.gtf).map{ it -> [ [id:'gtf'], it ] }.collect()
                                                : PREPARE_GENOME.out.gtf
    ch_star_index       = params.star_index     ? Channel.fromPath(params.star_index).map{ it -> [ [id:'star_index'], it ] }.collect()
                                                : PREPARE_GENOME.out.star_index
    ch_dbsnp_tbi        = params.dbsnp          ? params.dbsnp_tbi ? Channel.fromPath(params.dbsnp_tbi).collect()
                                                : PREPARE_GENOME.out.dbsnp_tbi : Channel.value([])
    ch_known_indels_tbi = params.known_indels   ? params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi).collect()
                                                : PREPARE_GENOME.out.known_indels_tbi : Channel.value([])

    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, rnavar will download the cache
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        snpeff_info     = Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    } else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_genome,
            params.snpeff_db,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            "Please refer to https://nf-co.re/rnavar/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

            snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
            vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    //
    // WORKFLOW: Run pipeline
    //
    RNAVAR (
        samplesheet
    )

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
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
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
