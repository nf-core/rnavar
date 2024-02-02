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
params.fasta             = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai         = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.dict              = WorkflowMain.getGenomeAttribute(params, 'dict')
params.gtf               = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff               = WorkflowMain.getGenomeAttribute(params, 'gff')
params.exon_bed          = WorkflowMain.getGenomeAttribute(params, 'exon_bed')
params.star_index        = WorkflowMain.getGenomeAttribute(params, 'star')
params.dbsnp             = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi         = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.known_indels      = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi  = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.snpeff_db         = WorkflowMain.getGenomeAttribute(params, 'snpeff_db')
params.snpeff_genome     = WorkflowMain.getGenomeAttribute(params, 'snpeff_genome')
params.vep_cache_version = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome        = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_species       = WorkflowMain.getGenomeAttribute(params, 'vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log, args)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNAVAR } from './workflows/rnavar'

//
// WORKFLOW: Run main nf-core/rnavar analysis pipeline
//
workflow NFCORE_RNAVAR {
    RNAVAR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_RNAVAR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
