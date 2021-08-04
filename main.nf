#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/rnavar
========================================================================================
    Github : https://github.com/nf-core/rnavar
    Website: https://nf-co.re/rnavar
    Slack  : https://nfcore.slack.com/channels/rnavar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta      = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf        = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff        = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed   = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.star_index = WorkflowMain.getGenomeAttribute(params, 'star')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { RNAVAR } from './workflows/rnavar'

//
// WORKFLOW: Run main nf-core/rnavar analysis pipeline
//
workflow NFCORE_RNAVAR {
    RNAVAR ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_RNAVAR ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
