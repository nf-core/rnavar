//
// Prepare reference genome files
//

include { BEDTOOLS_MERGE                    } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                     } from '../../../modules/nf-core/bedtools/sort/main'
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GFFREAD                           } from '../../../modules/nf-core/gffread/main'
include { GTF2BED                           } from '../../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate/main'
include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip/main'
include { TABIX_TABIX as TABIX_DBSNP        } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS } from '../../../modules/nf-core/tabix/tabix/main'

workflow PREPARE_GENOME {
    take:
    ch_fasta_raw      // file: /path/to/genome.fasta
    ch_gff            // file: /path/to/genome.gff
    ch_gtf_raw        // file: /path/to/genome.gtf
    ch_dbsnp
    ch_known_indels
    feature_type

    main:
    ch_versions = Channel.empty()

    //Unzip reference genome files if needed

    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA(ch_fasta_raw)

        ch_fasta = GUNZIP_FASTA.out.gunzip

    } else {
        ch_fasta = ch_fasta_raw
    }

    if (params.gtf.endsWith('.gz')) {
        GUNZIP_GTF(ch_gtf_raw)

        ch_gtf = GUNZIP_GTF.out.gunzip

    } else {
        ch_gtf = ch_gtf_raw
    }

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    GFFREAD(ch_gff, ch_fasta)
    SAMTOOLS_FAIDX(ch_fasta, [['id':'genome'], []])
    TABIX_DBSNP(ch_dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_INDELS(ch_known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )

    ch_gtf = ch_gtf.mix(GFFREAD.out.gtf)

    GTF2BED(ch_gtf, feature_type)
    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(GFFREAD.out.versions)
    ch_versions = ch_versions.mix(GTF2BED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)

    emit:
    dict             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              //    path: genome.fasta.dict
    exon_bed         = GTF2BED.out.bed.map{ bed -> [ [ id:bed.baseName ], bed ] }.collect() //    path: exon.bed
    fasta            = ch_fasta
    fasta_fai        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                     //    path: genome.fasta.fai
    gtf              = ch_gtf                                                               //    path: genome.gtf
    star_index       = STAR_GENOMEGENERATE.out.index                                        //    path: star/index/
    dbsnp_tbi        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()              // path: dbsnb.vcf.gz.tbi
    known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()       // path: {known_indels*}.vcf.gz.tbi
    versions         = ch_versions                                                          // channel: [ versions.yml ]
    // bedtools_sort    = ch_bedtools_sort    // path: sort.bed
    // bedtools_merge   = ch_bedtools_merge   // path: merge.bed
}
