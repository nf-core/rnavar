//
// This file holds several functions specific to the workflow/rnavar.nf in the nf-core/rnavar pipeline
//

class WorkflowRnavar {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.gtf && !params.gff) {
            log.error "No GTF or GFF3 annotation specified! The pipeline requires at least one of these files."
            System.exit(1)
        }

        if ((!params.skip_baserecalibration) && (!params.dbsnp && !params.known_indels)) {
            log.error "Known variants VCF file or its index is missing!. At least --dbsnp (and its index) or --known_indels (and its index) is required."
            System.exit(1)
        }

        if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('vep')) && (!params.genome || !params.vep_genome || !params.vep_species || !params.vep_cache_version)) {
            log.error "Species name (using --vep_species), genome assembly (either --genome or --vep_genome) and cache version (--vep_cache_version) are required to run VEP variant annotation."
            System.exit(1)
        }

        if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff')) && (!params.genome || !params.snpeff_db)) {
            log.error "Either --genome or --snpeff_db is required to run snpEff variant annotation."
            System.exit(1)
        }

    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}
