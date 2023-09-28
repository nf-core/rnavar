//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF           } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF           } from '../../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_STAR_INDEX      } from '../../../modules/nf-core/untar/main'

workflow PREPARE_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    // gtf                  //      file: /path/to/genome.gtf
    // gff                  //      file: /path/to/genome.gff
    // exon_bed             //      file: /path/to/gene.bed
    // prepare_tool_indices
    // feature_type

    main:
    ch_versions = Channel.empty()

    ch_fasta = fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }

    //
    // Uncompress genome fasta file if required
    //
    // if (fasta.endsWith('.gz')) {
    //     ch_fasta    = GUNZIP_FASTA([[:], fasta]).gunzip.map{ meta, fasta -> fasta }
    //     ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    // } else {
    //     ch_fasta = Channel.value(file(fasta))
    // }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    // if (gtf) {
    //     if (gtf.endsWith('.gz')) {
    //         ch_gtf      = GUNZIP_GTF([[:], gtf]).gunzip.map{ meta, gtf -> gtf }
    //         ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    //     } else ch_gtf = Channel.value(file(gtf))
    // } else if (gff) {
    //     if (gff.endsWith('.gz')) {
    //         ch_gff      = GUNZIP_GFF([[:], gff]).gunzip.map{ meta, gff -> gff }
    //         ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
    //     } else ch_gff = Channel.value(file(gff))
    //     ch_gtf      = GFFREAD(ch_gff).gtf
    //     ch_versions = ch_versions.mix(GFFREAD.out.versions)
    // }

    //
    // Uncompress exon BED annotation file or create from GTF if required
    //
    // if (exon_bed) {
    //     if (exon_bed.endsWith('.gz')) {
    //         exonGENE_BED(
    //             Channel.fromPath(exon_bed).map{ it -> [[id:it[0].baseName], it] }
    //         )
    //         ch_exon_bed = GUNZIP_GENE_BED.out.gunzip.map{ meta, bed -> [bed] }.collect()
    //         ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
    //     } else {
    //         ch_exon_bed = Channel.fromPath(exon_bed).collect()
    //     }
    // } else {
    //     ch_exon_bed = GTF2BED( ch_gtf , feature_type).bed.collect()
    //     ch_versions = ch_versions.mix(GTF2BED.out.versions)
    // }

    //ch_exon_bed.view()
    //ch_exon_bed.map{ it -> [[id:'exome'], it] }
    //ch_exon_bed.view()
    // Bedtools sort
    // ch_bedtools_sort = BEDTOOLS_SORT(ch_exon_bed.map{ it -> [[id:'exome'], it] }, 'sorted').sorted.collect()
    // ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)


    // // Bedtools merge
    // ch_bedtools_merge = BEDTOOLS_MERGE(ch_bedtools_sort).bed
    // ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)


    // Index the genome fasta
    SAMTOOLS_FAIDX(ch_fasta, [['id':null], []])

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Create dictionary file for the genome fasta
    // ch_fasta_dict = Channel.empty()
    // if (params.dict) ch_fasta_dict = Channel.fromPath(params.dict).collect()
    // else ch_fasta_dict = GATK4_CREATESEQUENCEDICTIONARY(ch_fasta).dict

    //
    // Uncompress STAR index or generate from scratch if required
    //
    // ch_star_index = Channel.empty()
    // if ('star' in prepare_tool_indices) {
    //     if (params.star_index) {
    //         if (params.star_index.endsWith('.tar.gz')) {
    //             UNTAR_STAR_INDEX(
    //                 Channel.fromPath(params.star_index).map{ it -> [[id:it[0].baseName], it] }
    //             )
    //             ch_star_index = UNTAR_STAR_INDEX.out.untar.map{ meta, star_index -> [star_index] }.collect()
    //             ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
    //         } else {
    //             ch_star_index = Channel.fromPath(params.star_index).collect()
    //         }
    //     }
    //     else {
    //         STAR_GENOMEGENERATE(
    //             ch_fasta,ch_gtf
    //         )
    //         .index
    //         .set { ch_star_index }
    //         ch_versions     = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    //     }
    // }


    emit:
    // fasta            = ch_fasta            // path: genome.fasta
    fasta_fai        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] } // path: genome.fasta.fai    dict             = ch_fasta_dict       // path: genome.fasta.dict
    // gtf              = ch_gtf              // path: genome.gtf
    // exon_bed         = ch_exon_bed         // path: exon.bed
    // bedtools_sort    = ch_bedtools_sort    // path: sort.bed
    // bedtools_merge   = ch_bedtools_merge   // path: merge.bed
    // star_index       = ch_star_index       // path: star/index/
    versions         = ch_versions         // channel: [ versions.yml ]
}
