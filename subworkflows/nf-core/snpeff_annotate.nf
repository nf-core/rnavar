//
// Run snpEff to annotate VCF files
//

//params.bgziptabix_snpeff = [:]
//params.snpeff_options    = [:]
//params.snpeff_tag        = [:]
//params.use_cache         = [:]

include { SNPEFF } from '../../modules/nf-core/modules/snpeff/main'
//addParams(
//    options:    params.snpeff_options,
//    snpeff_tag: params.snpeff_tag,
//    use_cache:  params.use_cache
//)

include { TABIX_BGZIPTABIX } from '../../modules/nf-core/modules/tabix/bgziptabix/main' //addParams(options: params.bgziptabix_snpeff_options)

workflow SNPEFF_ANNOTATE {
    take:
    vcf            // channel: [ val(meta), vcf, tbi ]
    snpeff_db      // value: version of db to use
    snpeff_cache   // path: path_to_snpeff_cache (optionnal)

    main:

    ch_versions = Channel.empty()

    SNPEFF(vcf, snpeff_db, snpeff_cache)
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    vcf            = TABIX_BGZIPTABIX.out.tbi // channel: [ val(meta), vcf, tbi ]
    snpeff_report  = SNPEFF.out.report        // path: *.html
    versions       = ch_versions              // channel: [versions.yml]
}
