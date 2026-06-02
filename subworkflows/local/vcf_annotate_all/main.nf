//
// ANNOTATION
//

include { BCFTOOLS_ANNOTATE                             } from '../../../modules/nf-core/bcftools/annotate'
include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_SNPEFF                           } from '../../nf-core/vcf_annotate_snpeff'

workflow VCF_ANNOTATE_ALL {
    take:
    vcf // channel: [ val(meta), vcf ] (mandatory)
    fasta // channel: [ val(meta), fasta ] (optional)
    tools // array: list of tools to apply (mandatory)
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files
    bcftools_annotations
    bcftools_annotations_index
    bcftools_columns
    bcftools_header_lines

    main:
    json_ann = channel.empty()
    tab_ann = channel.empty()
    vcf_ann = channel.empty()

    if (bcftools_columns) {
        vcf_for_bcfann = vcf
            .combine(bcftools_annotations)
            .combine(bcftools_annotations_index)
            .combine(bcftools_columns)
            .combine(bcftools_header_lines)
            .map { meta, vcf_, annotation, annotation_index, columns, header_file -> [meta, vcf_, [], annotation, annotation_index, columns, header_file, []] }
    }
    else {
        vcf_for_bcfann = vcf
            .combine(bcftools_annotations)
            .combine(bcftools_annotations_index)
            .combine(bcftools_header_lines)
            .map { meta, vcf_, annotation, annotation_index, header_file -> [meta, vcf_, [], annotation, annotation_index, [], header_file, []] }
    }

    BCFTOOLS_ANNOTATE(vcf_for_bcfann.filter { 'bcfann' in tools })
    vcf_ann = vcf_ann.mix(BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.index, failOnDuplicate: true, failOnMismatch: true))

    VCF_ANNOTATE_SNPEFF(vcf.filter { ('merge' in tools || 'snpeff' in tools) }, snpeff_db, snpeff_cache)
    vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)

    VCF_ANNOTATE_MERGE(
        VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map { meta, vcf_, _tbi -> [meta, vcf_, []] }.filter { 'merge' in tools },
        fasta,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        vep_extra_files,
    )

    vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
    tab_ann = tab_ann.mix(VCF_ANNOTATE_MERGE.out.tab)
    json_ann = json_ann.mix(VCF_ANNOTATE_MERGE.out.json)

    VCF_ANNOTATE_ENSEMBLVEP(vcf.map { meta, vcf_ -> [meta, vcf_, []] }.filter { 'vep' in tools }, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

    vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
    tab_ann = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
    json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)

    emit:
    vcf_ann // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
}
