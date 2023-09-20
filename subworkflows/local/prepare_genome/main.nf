//
// Prepare reference genome files
//

include { BEDTOOLS_MERGE                 } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                  } from '../../../modules/nf-core/bedtools/sort/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GFFREAD                        } from '../../../modules/nf-core/gffread/main'
include { GTF2BED                        } from '../../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE            } from '../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME {
    take:
    ch_exon_bed   // file: /path/to/gene.bed
    ch_fasta      // file: /path/to/genome.fasta
    ch_gff        // file: /path/to/genome.gff
    ch_gtf        // file: /path/to/genome.gtf
    feature_type

    main:
    ch_versions = Channel.empty()

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    GFFREAD(ch_gff)
    SAMTOOLS_FAIDX(ch_fasta, [['id':null], []])

    ch_gtf = ch_gtf.mix(GFFREAD.out.gtf)

    GTF2BED(ch_gtf, feature_type)
    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(GFFREAD.out.versions)
    ch_versions = ch_versions.mix(GTF2BED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    dict       = GATK4_CREATESEQUENCEDICTIONARY.out.dict          //    path: genome.fasta.dict
    exon_bed   = GTF2BED.out.bed.collect()                        //    path: exon.bed
    fasta_fai  = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] } //    path: genome.fasta.fai
    gtf        = ch_gtf                                           //    path: genome.gtf
    star_index = STAR_GENOMEGENERATE.out.index                    //    path: star/index/
    versions   = ch_versions                                      // channel: [ versions.yml ]
    // bedtools_sort    = ch_bedtools_sort    // path: sort.bed
    // bedtools_merge   = ch_bedtools_merge   // path: merge.bed
}
