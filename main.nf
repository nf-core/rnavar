#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnavar
========================================================================================
 nf-core/rnavar Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnavar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta        = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf          = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff          = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed     = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.star_index   = WorkflowMain.getGenomeAttribute(params, 'star')


/*
========================================================================================
VALIDATE AND PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
RUN THE WORKFLOW
========================================================================================
*/

workflow {
    include { RNASEQ_VAR } from './workflows/rnavar'
    RNASEQ_VAR ()
  
}

/*
END OF THE FILE
*/
