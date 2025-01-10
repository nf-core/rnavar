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
    fasta_raw      // file: /path/to/genome.fasta
    gff            // file: /path/to/genome.gff
    gtf_raw        // file: /path/to/genome.gtf
    dbsnp
    known_indels
    feature_type

    main:
    ch_versions = Channel.empty()

    //Unzip reference genome files if needed

    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA(fasta_raw)

        fasta = GUNZIP_FASTA.out.gunzip

    } else {
        fasta = fasta_raw
    }

    if (params.gtf.endsWith('.gz')) {
        GUNZIP_GTF(gtf_raw)

        gtf = GUNZIP_GTF.out.gunzip

    } else {
        gtf = gtf_raw
    }

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    GFFREAD(gff, fasta)
    SAMTOOLS_FAIDX(fasta, [['id':'genome'], []])
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )

    gtf = gtf.mix(GFFREAD.out.gtf)

    GTF2BED(gtf, feature_type)
    STAR_GENOMEGENERATE(fasta, gtf)

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
    fasta            = fasta
    fasta_fai        = SAMTOOLS_FAIDX.out.fai                                               //    path: genome.fasta.fai
    gtf              = gtf.first()                                                       //    path: genome.gtf
    star_index       = STAR_GENOMEGENERATE.out.index.first()                                //    path: star/index/
    dbsnp_tbi        = TABIX_DBSNP.out.tbi.map{ _meta, tbi -> [tbi] }.collect()              // path: dbsnb.vcf.gz.tbi
    known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map{ _meta, tbi -> [tbi] }.collect()       // path: {known_indels*}.vcf.gz.tbi
    versions         = ch_versions                                                          // channel: [ versions.yml ]
    // bedtools_sort    = bedtools_sort    // path: sort.bed
    // bedtools_merge   = bedtools_merge   // path: merge.bed
}
