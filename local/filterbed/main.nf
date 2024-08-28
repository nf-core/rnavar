process FILTERBEDFILE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

     input:
    tuple val(meta), path(bed)
    path  dict_file

    output:
    tuple val(meta), path('filtered.bed'), emit: filtered_bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<PYCODE
    import argparse
    import os

    def load_sequences_from_dict(dict_file):
        sequences = set()
        with open(dict_file, 'r') as file:
            for line in file:
                if line.startswith('@SQ'):
                    parts = line.split('\t')
                    for part in parts:
                        if part.startswith('SN:'):
                            sequences.add(part.split(':')[1])
        return sequences

    def filter_bed_file(bed_file, sequences, output_file):
        with open(bed_file, 'r') as file, open(output_file, 'w') as out:
            for line in file:
                sequence = line.split('\t')[0]
                if sequence in sequences:
                    out.write(line)

    def main(bed_file, dict_file, output_file):
        sequences = load_sequences_from_dict(dict_file)
        filter_bed_file(bed_file, sequences, output_file)
        print(f"Output file {output_file} created in {os.getcwd()}")

    if __name__ == "__main__":
        main("${bed}", "${dict_file}", "filtered.bed")

    PYCODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
