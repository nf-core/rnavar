// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]

process GTF2BED {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bedops=2.4.39" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedops%3A2.4.39--h7d875b9_1"
    } else {
        container "quay.io/biocontainers/bedops:2.4.39--h7d875b9_1"
    }

    input:
    path gtf

    output:
    path "*.bed"       , emit: bed
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def grep_interval_featureType  = params.interval_list_featureType ? "| grep -iw '$params.interval_list_featureType' | cut -f1-3" : ""
    """
    gtf2bed < \\
        $gtf $grep_interval_featureType \\
        > ${gtf.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        bedops: \$(echo \$(bedops --version 2>&1) | sed 's/.*version: \\(.*\\) (typical).*/\\1/')
    END_VERSIONS
    """
}
