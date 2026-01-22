/**
 * Aggregate results from multiple analysis tools into a single report.
 *
 * MultiQC searches a given directory for analysis logs and compiles them
 * into a single HTML report. It supports output from many common
 * bioinformatics tools including FastQC, STAR, Picard, GATK, and more.
 *
 * The report provides:
 * - Summary statistics across all samples
 * - Interactive plots for QC metrics
 * - Data tables for detailed metrics
 * - Export functionality for plots and data
 *
 * @param multiqc_files Input files/directories to search for tool outputs
 * @param multiqc_config Custom MultiQC configuration file (optional)
 * @param extra_multiqc_config Additional MultiQC config to merge (optional)
 * @param multiqc_logo Custom logo for the report header (optional)
 * @param replace_names File with sample name replacements (optional)
 * @param sample_names File with sample name mappings (optional)
 */
process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data' :
        'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b' }"

    input:
    path multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "*.html"      , emit: report
    path "*_data"      , emit: data
    path "*_plots"     , optional:true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    def extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${extra_config} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
