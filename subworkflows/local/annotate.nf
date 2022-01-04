//
// ANNOTATION
//

include { SNPEFF_ANNOTATE } from '../nf-core/snpeff_annotate'
include { ENSEMBLVEP_ANNOTATE as MERGE_ANNOTATE } from '../nf-core/ensemblvep_annotate'
include { ENSEMBLVEP_ANNOTATE } from '../nf-core/ensemblvep_annotate'

workflow ANNOTATE {
    take:
    vcf          // channel: [ val(meta), vcf ]
    tools
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache

    main:

    ch_versions       = Channel.empty()
    merge_vcf_ann     = Channel.empty()
    merge_vep_report  = Channel.empty()
    merge_vep_version = Channel.empty()
    snpeff_report     = Channel.empty()
    snpeff_vcf_ann    = Channel.empty()
    snpeff_version    = Channel.empty()
    vep_report        = Channel.empty()
    vep_vcf_ann       = Channel.empty()
    vep_version       = Channel.empty()

    if ('snpeff' in tools || 'merge' in tools) {
        (snpeff_vcf_ann, snpeff_report, snpeff_version) = SNPEFF_ANNOTATE(vcf, snpeff_db, snpeff_cache)
        ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions.first())
    }

    if ('merge' in tools) {
        vcf_ann_for_merge = snpeff_vcf_ann.map{ meta, vcf, tbi -> [meta, vcf, tbi] }
        (merge_vcf_ann, merge_vep_report, merge_vep_version) = MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE.out.versions.first())
    }

    if ('vep' in tools) {
        (vep_vcf_ann, vep_report, vep_version) = ENSEMBLVEP_ANNOTATE(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
        ch_versions = ch_versions.mix(ENSEMBLVEP_ANNOTATE.out.versions.first())
    }

    vcf_ann     = snpeff_vcf_ann.mix(merge_vcf_ann, vep_vcf_ann)
    reports     = snpeff_report.mix(merge_vep_report, vep_report)
    versions    = ch_versions

    emit:
        reports
        vcf_ann
        versions
}
