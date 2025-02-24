//
// Prepare reference genome files
//

include { BEDTOOLS_MERGE                                } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                                 } from '../../../modules/nf-core/bedtools/sort/main'
include { GATK4_CREATESEQUENCEDICTIONARY                } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GFFREAD                                       } from '../../../modules/nf-core/gffread/main'
include { GTF2BED                                       } from '../../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX                                } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE                           } from '../../../modules/nf-core/star/genomegenerate/main'
include { GUNZIP as GUNZIP_FASTA                        } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                          } from '../../../modules/nf-core/gunzip/main'
include { TABIX_TABIX as TABIX_DBSNP                    } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIPTABIX_DBSNP          } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS             } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIPTABIX_KNOWN_INDELS   } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow PREPARE_GENOME {
    take:
    fasta_raw           // channel: /path/to/genome.fasta
    gff                 // channel: /path/to/genome.gff
    gtf_raw             // channel: /path/to/genome.gtf
    dbsnp               // channel: /path/to/dbnsp.vcf.gz
    known_indels        // channel: [/path/to/known_indels]
    known_indels_tbi    // channel: [/path/to/known_indels_index]
    feature_type

    main:
    def ch_versions = Channel.empty()

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

    def ch_dbsnp = Channel.value([[],[]])
    def ch_dbsnp_tbi = Channel.value([[],[]])
    if (params.dbsnp && !params.dbsnp.toString().endsWith(".gz")) {
        BGZIPTABIX_DBSNP(
            dbsnp
        )
        ch_versions = ch_versions.mix(BGZIPTABIX_DBSNP.out.versions.first())
        ch_dbsnp = BGZIPTABIX_DBSNP.out.gz_tbi.map { meta, vcf, _tbi -> [ meta, vcf ] }.collect()
        ch_dbsnp_tbi = BGZIPTABIX_DBSNP.out.gz_tbi.map { meta, _vcf, tbi -> [ meta, tbi ] }.collect()
    } else if(params.dbsnp && !params.dbsnp_tbi) {
        TABIX_DBSNP(
            dbsnp
        )
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions.first())
        ch_dbsnp_tbi = TABIX_DBSNP.out.tbi.collect()
    }

    def ch_known_indels_input = known_indels.map { file -> [[id:file.name], file ] }
        .join(known_indels_tbi.map { file -> [[id:file.baseName], file ] }, failOnDuplicate:true, remainder:true)
        .branch { meta, file, index ->
            plain: !file.toString().endsWith(".gz")
                return [ meta, file ]
            bgzip_noindex: !index && file.toString().endsWith(".gz")
                return [ meta, file ]
            bgzip_index: true
                return [ meta, file, index ]
        }

    BGZIPTABIX_KNOWN_INDELS(
        ch_known_indels_input.plain
    )
    ch_versions = ch_versions.mix(BGZIPTABIX_KNOWN_INDELS.out.versions.first())

    TABIX_KNOWN_INDELS(
        ch_known_indels_input.bgzip_noindex
    )
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions.first())

    def ch_known_indels = BGZIPTABIX_KNOWN_INDELS.out.gz_tbi.map { _meta, file, _index -> file }
        .mix(ch_known_indels_input.bgzip_noindex.map { _meta, file -> file })
        .mix(ch_known_indels_input.bgzip_index.map { _meta, file, _index -> file })
        .collect()

    def ch_known_indels_tbi = BGZIPTABIX_KNOWN_INDELS.out.gz_tbi.map { _meta, _file, index -> index }
        .mix(TABIX_KNOWN_INDELS.out.tbi.map { _meta, tbi -> tbi })
        .mix(ch_known_indels_input.bgzip_index.map { _meta, _file, index -> index })
        .collect()

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    GFFREAD(gff, fasta)
    SAMTOOLS_FAIDX(fasta, [['id':'genome'], []])

    gtf = gtf.mix(GFFREAD.out.gtf)

    GTF2BED(gtf, feature_type)
    STAR_GENOMEGENERATE(fasta, gtf)

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(GFFREAD.out.versions)
    ch_versions = ch_versions.mix(GTF2BED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    dict             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              //    path: genome.fasta.dict
    exon_bed         = GTF2BED.out.bed.map{ bed -> [ [ id:bed.baseName ], bed ] }.collect() //    path: exon.bed
    fasta            = fasta
    fasta_fai        = SAMTOOLS_FAIDX.out.fai                                               //    path: genome.fasta.fai
    gtf              = gtf.first()                                                       //    path: genome.gtf
    star_index       = STAR_GENOMEGENERATE.out.index.first()                                //    path: star/index/
    dbsnp            = ch_dbsnp                                                               // path: dbsnp.vcf.gz
    dbsnp_tbi        = ch_dbsnp_tbi                                                           // path: dbsnp.vcf.gz.tbi
    known_indels     = ch_known_indels                                                      // path: {known_indels*}.vcf.gz
    known_indels_tbi = ch_known_indels_tbi                                                 // path: {known_indels*}.vcf.gz.tbi
    versions         = ch_versions                                                          // channel: [ versions.yml ]
    // bedtools_sort    = bedtools_sort    // path: sort.bed
    // bedtools_merge   = bedtools_merge   // path: merge.bed
}
