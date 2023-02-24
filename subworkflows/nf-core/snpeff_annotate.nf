//
// Run snpEff to annotate VCF files
//

include { SNPEFF_SNPEFF } from '../../modules/nf-core/snpeff/snpeff/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow SNPEFF_ANNOTATE {
    take:
    vcf            // channel: [ val(meta), vcf, tbi ]
    snpeff_db      // value: version of db to use
    snpeff_cache   // path: path_to_snpeff_cache (optionnal)

    main:

    ch_versions = Channel.empty()

    SNPEFF_SNPEFF (
        vcf,
        snpeff_db,
        snpeff_cache
    )
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions.first())

    TABIX_BGZIPTABIX (
        SNPEFF_SNPEFF.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi     = TABIX_BGZIPTABIX.out.gz_tbi    // channel: [ val(meta), vcf, tbi ]
    reports     = SNPEFF_SNPEFF.out.report              // path: *.html
    versions    = ch_versions                    // channel: [versions.yml]
}
