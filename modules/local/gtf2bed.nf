process GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "conda-forge::r-base=3.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:3.5.0' :
        'biocontainers/r-base:3.5.0'}"

    input:
    path gtf
    val feature_type
    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def allowed_type = ["exon", "transcript", "gene"];
    if (feature_type){
        feature_type = allowed_type.contains(feature_type) ? feature_type : "exon"
    }
    """
    Rscript --no-save -<<'RCODE'
        gtf = read.table("${gtf}", sep="\t")
        gtf = subset(gtf, V3 == "${feature_type}")
        write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "tmp.exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
    RCODE

    awk '{print \$1 "\t" (\$2 - 1) "\t" \$3}' tmp.exome.bed > exome.bed
    rm -rf tmp.exome.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(echo \$(Rscript --version 2>&1) | sed 's/R scripting front-end version //')
    END_VERSIONS
    """
}
