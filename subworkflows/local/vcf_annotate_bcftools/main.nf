
//
// Run BCFtools to annotate VCF files
//

include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_BCFTOOLS {
    take:
    ch_vcf               // channel: [ val(meta), vcf ]

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_ANNOTATE(ch_vcf.map { meta, vcf -> [ meta, vcf, [], [], [], []]})
    TABIX_TABIX(BCFTOOLS_ANNOTATE.out.vcf)

    ch_vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    versions = ch_versions //    path: versions.yml
}
