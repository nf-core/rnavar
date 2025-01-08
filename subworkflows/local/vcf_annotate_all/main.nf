//
// ANNOTATION
//

include { VCF_ANNOTATE_BCFTOOLS                         } from '../vcf_annotate_bcftools/main'
include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { VCF_ANNOTATE_SNPEFF                           } from '../../nf-core/vcf_annotate_snpeff/main'

workflow VCF_ANNOTATE_ALL {
    take:
    ch_vcf          // channel: [ val(meta), vcf ]
    ch_fasta
    tools        // Mandatory, list of tools to apply
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files

    main:
    def ch_reports = Channel.empty()
    def vcf_ann = Channel.empty()
    def tab_ann = Channel.empty()
    def json_ann = Channel.empty()
    def versions = Channel.empty()

    if (tools.split(',').contains('bcfann')) {
        VCF_ANNOTATE_BCFTOOLS(ch_vcf)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_BCFTOOLS.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_BCFTOOLS.out.versions)
    }


    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        VCF_ANNOTATE_SNPEFF(ch_vcf, snpeff_db, snpeff_cache)

        ch_reports = ch_reports.mix(VCF_ANNOTATE_SNPEFF.out.reports.map{ _meta, reports -> [ reports ] })
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    }

    if (tools.split(',').contains('merge')) {
        vcf_ann_for_merge = VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map{ meta, vcf, _tbi -> [ meta, vcf, [] ] }
        VCF_ANNOTATE_MERGE(vcf_ann_for_merge, ch_fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports = ch_reports.mix(VCF_ANNOTATE_MERGE.out.reports)
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    }

    if (tools.split(',').contains('vep')) {
        vcf_for_vep = ch_vcf.map{ meta, vcf -> [ meta, vcf, [] ] }
        VCF_ANNOTATE_ENSEMBLVEP(vcf_for_vep, ch_fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports  = ch_reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_ann  = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
        versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    }

    emit:
    vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports = ch_reports      //    path: *.html
    versions     //    path: versions.yml
}
