//
// Uncompress and prepare reference genome files
//

include { BEDTOOLS_MERGE                    } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                     } from '../../../modules/nf-core/bedtools/sort/main'
include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GFFREAD                           } from '../../../modules/nf-core/gffread/main'
include { GTF2BED                           } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                        } from '../../../modules/local/gtf_filter'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate/main'
include { TABIX_TABIX as TABIX_DBSNP        } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    gene_bed                 //      file: /path/to/gene.bed
    dbsnp                    //      file: /path/to/dbsnp.vcf
    known_indels             //      file: /path/to/known_indels.vcf
    star_index               // directory: /path/to/star/index/
    featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    gencode                  //   boolean: whether the genome is from GENCODE
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs

    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff))
            }
            ch_gtf      = GFFREAD ( ch_gff ).gtf
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        // Determine whether to filter the GTF or not
        def filter_gtf =
            ((
                // Condition 1: Transcript FASTA file is not provided
                !transcript_fasta
            )) &&
            (
                // Condition 2: --skip_gtf_filter is not provided
                !skip_gtf_filter
            )
        if (filter_gtf) {
            GTF_FILTER ( ch_fasta, ch_gtf )
            ch_gtf = GTF_FILTER.out.genome_gtf
            ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        }
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : featurecounts_group_type
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta))
        }

        CUSTOM_CATADDITIONALFASTA (
            ch_fasta.combine(ch_gtf).map { fasta, gtf -> [ [:], fasta, gtf ] },
            ch_add_fasta.map { [ [:], it ] },
            biotype
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map { it[1] }.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map { it[1] }.first()
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed))
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    TABIX_DBSNP(ch_dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_INDELS(ch_known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (star_index) {
        if (star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
            ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = Channel.value(file(star_index))
        }
    } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    emit:
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    dbsnp_tbi        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()              // path: dbsnb.vcf.gz.tbi
    dict             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              //    path: genome.fasta.dict
    fai              = ch_fai                    // channel: path(genome.fai)
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()       // path: {known_indels*}.vcf.gz.tbi
    star_index       = ch_star_index             // channel: path(star/index/)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
