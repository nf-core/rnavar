
//
// Run BCFtools to annotate VCF files
//

include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_BCFTOOLS {
    take:
    vcf               // channel: [ val(meta), vcf ]

    main:
    def ch_versions = Channel.empty()

    BCFTOOLS_ANNOTATE(vcf.map { meta, vcf_ -> [ meta, vcf_, [], [], [], []]})
    TABIX_TABIX(BCFTOOLS_ANNOTATE.out.vcf)

    vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf.join(TABIX_TABIX.out.tbi, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = vcf_tbi  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    versions = ch_versions //    path: versions.yml
}
