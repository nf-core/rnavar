//
// ANNOTATION
//

include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { VCF_ANNOTATE_SNPEFF                           } from '../../nf-core/vcf_annotate_snpeff/main'

workflow VCF_ANNOTATE_ALL {
    take:
    vcf          // channel: [ val(meta), vcf ]
    fasta
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
    def vcf_annotated = Channel.empty()
    def tab_annotated = Channel.empty()
    def json_annotated = Channel.empty()
    def ch_versions = Channel.empty()

    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        VCF_ANNOTATE_SNPEFF(vcf, snpeff_db, snpeff_cache)

        ch_reports = ch_reports.mix(VCF_ANNOTATE_SNPEFF.out.reports.map{ _meta, reports_ -> [ reports_ ] })
        vcf_annotated = vcf_annotated.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
        ch_versions = ch_versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    }

    if (tools.split(',').contains('merge')) {
        vcf_ann_for_merge = VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map{ meta, vcf_, _tbi -> [ meta, vcf_, [] ] }
        VCF_ANNOTATE_MERGE(vcf_ann_for_merge, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports = ch_reports.mix(VCF_ANNOTATE_MERGE.out.reports)
        vcf_annotated = vcf_annotated.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        ch_versions = ch_versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    }

    if (tools.split(',').contains('vep')) {
        vcf_for_vep = vcf.map{ meta, vcf_ -> [ meta, vcf_, [] ] }
        VCF_ANNOTATE_ENSEMBLVEP(vcf_for_vep, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports  = ch_reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
        vcf_annotated  = vcf_annotated.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_annotated  = tab_annotated.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_annotated = json_annotated.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
        ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    }

    emit:
    vcf_ann = vcf_annotated     // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann = tab_annotated
    json_ann = json_annotated
    reports = ch_reports      //    path: *.html
    versions = ch_versions     //    path: ch_versions.yml
}
