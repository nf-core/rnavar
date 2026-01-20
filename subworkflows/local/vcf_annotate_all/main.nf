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

    if (tools.split(',').contains('bcfann')) {
        BCFTOOLS_ANNOTATE(
            vcf.map { meta, vcf_ -> [meta, vcf_, []] }.combine(bcftools_annotations).combine(bcftools_annotations_index),
            bcftools_columns,
            bcftools_header_lines,
            [],
        )

        vcf_ann = vcf_ann.mix(BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.tbi, failOnDuplicate: true, failOnMismatch: true))
    }


    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        VCF_ANNOTATE_SNPEFF(vcf, snpeff_db, snpeff_cache)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
    }

    if (tools.split(',').contains('merge')) {
        vcf_ann_for_merge = VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map { meta, vcf_, _tbi -> [meta, vcf_, []] }
        VCF_ANNOTATE_MERGE(vcf_ann_for_merge, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        tab_ann = tab_ann.mix(VCF_ANNOTATE_MERGE.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_MERGE.out.json)
    }

    if (tools.split(',').contains('vep')) {
        vcf_for_vep = vcf.map { meta, vcf_ -> [meta, vcf_, []] }
        VCF_ANNOTATE_ENSEMBLVEP(vcf_for_vep, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_ann = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
    }

    emit:
    vcf_ann // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
}
