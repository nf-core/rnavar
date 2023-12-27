process FILTERBEDFILE {
    tag "${genome_bed} -> ${filtered_bed}"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path genome_bed
    path genome_dict

    output:
    path 'filtered_exome.bed', emit: filtered_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<PYCODE
    import argparse

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

    def main():
        sequences = load_sequences_from_dict("${genome_dict}")
        filter_bed_file("${genome_bed}", sequences, "filtered_exome.bed")

    if __name__ == "__main__":
        main()
    PYCODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
