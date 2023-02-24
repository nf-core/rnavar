//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow ENSEMBLVEP_ANNOTATE {
    take:
    vcf               // channel: [ val(meta), vcf, tbi ]
    vep_genome        //   value: which genome
    vep_species       //   value: which species
    vep_cache_version //   value: which cache version
    vep_cache         //    path: path_to_vep_cache (optionnal)

    main:

    ch_versions = Channel.empty()

    ENSEMBLVEP_VEP (
        vcf,
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())

    TABIX_BGZIPTABIX (
        ENSEMBLVEP_VEP.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi     = TABIX_BGZIPTABIX.out.gz_tbi    // channel: [ val(meta), vcf, tbi ]
    reports     = ENSEMBLVEP_VEP.out.report          // path: *.html
    versions    = ch_versions                    // channel: [versions.yml]
}
