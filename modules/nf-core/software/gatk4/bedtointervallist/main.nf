nextflow.enable.dsl=2

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_BEDTOINTERVALLIST {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::gatk4=4.1.9.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0'
    } else {
        container 'quay.io/biocontainers/gatk4:4.1.9.0--py39_0'
    }

    input:
    path bed
    path sequence_dict

    output:
    path '*.interval_list'                 , emit: interval_list
    path '*.version.txt'                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${bed}.${options.suffix}" : "${bed}"
    """
    gatk BedToIntervalList \\
        -I $bed \\
        -SD $sequence_dict \\
        -O ${prefix}.interval_list \\
        $options.args

    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
