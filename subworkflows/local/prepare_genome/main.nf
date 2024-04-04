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
include { UNZIP as UNZIP_FASTA           } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GTF             } from '../../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GFF             } from '../../../modules/nf-core/unzip/main'

workflow PREPARE_GENOME {
    take:
    ch_fasta_raw      // file: /path/to/genome.fasta
    ch_gff            // file: /path/to/genome.gff
    ch_gtf_raw        // file: /path/to/genome.gtf
    feature_type

    main:
    ch_versions = Channel.empty()

    //Unzip reference genome files if needed

    if (params.fasta.endsWith('.gz')) {
        UNZIP_FASTA(ch_fasta_raw)

        //file gets saved into a folder with the same name as the file, need to add the missing depth to the folder
        ch_fasta = UNZIP_FASTA.out.unzipped_archive.map{ meta, file ->
            def file_name      = file.baseName
            def full_file_path = file.toString()  + '/' + file_name + '.fa'

            [meta, full_file_path] }
        } else {
        ch_fasta = ch_fasta_raw
    }

    if (params.gtf.endsWith('.gz')) {
        UNZIP_GTF(ch_gtf_raw)

        //file gets saved into a folder with the same name as the file, need to add the missing depth to the folder
        ch_gtf = UNZIP_GTF.out.unzipped_archive.map{ meta, file ->
            def file_name      = file.baseName
            def full_file_path = file.toString()  + '/' + file_name + '.gtf'

            [meta,full_file_path] }
        } else {
        ch_gtf = ch_gtf_raw
    }

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    GFFREAD(ch_gff)
    SAMTOOLS_FAIDX(ch_fasta, [['id':'genome'], []])

    ch_gtf = ch_gtf.mix(GFFREAD.out.gtf)

    GTF2BED(ch_gtf, feature_type)
    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(GFFREAD.out.versions)
    ch_versions = ch_versions.mix(GTF2BED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    dict       = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              //    path: genome.fasta.dict
    exon_bed   = GTF2BED.out.bed.map{ bed -> [ [ id:bed.baseName ], bed ] }.collect() //    path: exon.bed
    fasta      = ch_fasta
    fasta_fai  = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                     //    path: genome.fasta.fai
    gtf        = ch_gtf                                                               //    path: genome.gtf
    star_index = STAR_GENOMEGENERATE.out.index                                        //    path: star/index/
    versions   = ch_versions                                                          // channel: [ versions.yml ]
    // bedtools_sort    = ch_bedtools_sort    // path: sort.bed
    // bedtools_merge   = ch_bedtools_merge   // path: merge.bed
}
