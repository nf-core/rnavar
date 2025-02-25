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
include { UNTAR                                         } from '../../../modules/nf-core/untar'
include { STAR_INDEXVERSION                             } from '../../../modules/nf-core/star/indexversion'

workflow PREPARE_GENOME {
    take:
    fasta_raw           // channel: /path/to/genome.fasta
    star_index          // channel: /path/to/star_index
    gff                 // channel: /path/to/genome.gff
    gtf_raw             // channel: /path/to/genome.gtf
    dbsnp               // channel: /path/to/dbnsp.vcf.gz
    known_indels        // channel: [/path/to/known_indels]
    known_indels_tbi    // channel: [/path/to/known_indels_index]
    feature_type

    main:
    def ch_versions = Channel.empty()

    //Unzip reference genome files if needed

    def ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA(fasta_raw)

        ch_fasta = GUNZIP_FASTA.out.gunzip

    } else {
        ch_fasta = fasta_raw
    }

    def ch_gtf = Channel.empty()
    if (params.gtf.toString().endsWith('.gz')) {
        GUNZIP_GTF(gtf_raw)
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        ch_gtf = GUNZIP_GTF.out.gunzip
    } else if(params.gff) {
        GFFREAD(gff, ch_fasta.map { _meta, fasta_ -> fasta_}.collect())
        ch_versions = ch_versions.mix(GFFREAD.out.versions)

        ch_gtf = GFFREAD.out.gtf
    } else {
        ch_gtf = gtf_raw
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

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    SAMTOOLS_FAIDX(ch_fasta, [['id':'genome'], []])


    GTF2BED(ch_gtf, feature_type)

    //
    // STAR index handling
    //

    def star_index_input = star_index.branch { _meta, index ->
        tarzipped: index.name.endsWith(".tar.gz")
        no_index: !index
        index: true
    }

    UNTAR(
        star_index_input.tarzipped
    )
    ch_versions = ch_versions.mix(UNTAR.out.versions)

    STAR_INDEXVERSION()
    ch_versions = ch_versions.mix(STAR_INDEXVERSION.out.versions)

    def star_index_check = star_index_input.index
        .mix(UNTAR.out.untar)
        .combine(STAR_INDEXVERSION.out.index_version)
        .branch { meta, index, version_file ->
            def is_compatible = true
            if(!workflow.stubRun) {
                def minimal_version = version_file.text.replace("\n", "")
                def index_version = file("${index.toUri()}/genomeParameters.txt", checkIfExists:true)
                    .text
                    .readLines()
                    .find { line -> line.startsWith("versionGenome") }
                    .tokenize("\t")[-1]
                is_compatible = isCompatibleStarIndex(index_version, minimal_version)
                if(!is_compatible) {
                    log.warn("Detected a wrong version of the STAR index, expected a minimum version of ${minimal_version}. Automatically recreating the index of STAR...")
                }
            }
            compatible: is_compatible
                return [ meta, index]
            incompatible: !is_compatible
                return [ meta, [] ]
        }

    def genomegenerate_input = star_index_check.incompatible
        .mix(star_index_input.no_index)
        .combine(ch_fasta)
        .map { _meta1, _wrong_index, meta2, fasta ->
            [ meta2, fasta ]
        }

    STAR_GENOMEGENERATE(genomegenerate_input, ch_gtf)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    def star_index_output = STAR_GENOMEGENERATE.out.index
        .mix(star_index_check.compatible)

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(GTF2BED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    dict             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                              //    path: genome.fasta.dict
    exon_bed         = GTF2BED.out.bed.map{ bed -> [ [ id:bed.baseName ], bed ] }.collect() //    path: exon.bed
    fasta            = ch_fasta
    fasta_fai        = SAMTOOLS_FAIDX.out.fai                                               //    path: genome.fasta.fai
    gtf              = ch_gtf                                                                  //    path: genome.gtf
    star_index       = STAR_GENOMEGENERATE.out.index.first()                                //    path: star/index/
    dbsnp            = ch_dbsnp                                                               // path: dbsnp.vcf.gz
    dbsnp_tbi        = ch_dbsnp_tbi                                                           // path: dbsnp.vcf.gz.tbi
    known_indels     = ch_known_indels                                                      // path: {known_indels*}.vcf.gz
    known_indels_tbi = ch_known_indels_tbi                                                 // path: {known_indels*}.vcf.gz.tbi
    versions         = ch_versions                                                          // channel: [ versions.yml ]
    // bedtools_sort    = bedtools_sort    // path: sort.bed
    // bedtools_merge   = bedtools_merge   // path: merge.bed
}

def isCompatibleStarIndex(index_version, minimal_index_version) {
    def is_compatible = true
    if(minimal_index_version.isNumber()) {
        // Older version of STAR used a numerical versioning.
        // Return true if the index doesn't use the numerical versioning anymore
        if(!index_version.isNumber()) {
            is_compatible = true
        } else {
            is_compatible = index_version.toInteger() >= minimal_index_version.toInteger()
        }
    } else {
        if(index_version.isNumber()) {
            is_compatible = false
        } else {
            // Correctly compare semantic version strings: e.g 2.7.11b > 2.7.4a
            def min_list = convertVersionToList(minimal_index_version)
            def ind_list = convertVersionToList(index_version)
            ind_list.eachWithIndex { digit, idx ->
                if(digit > min_list[idx]) {
                    is_compatible = true
                    return
                } else if(digit < min_list[idx]) {
                    is_compatible = false
                    return
                }
            }
        }
    }
    return is_compatible
}

def convertVersionToList(version) {
    def init_list = version.tokenize(".")
    if(!init_list[-1].isNumber()) {
        // Handle cases where the last digit in the version contains a character: e.g. 2.7.11b
        def last_digit = init_list[-1]
        def numbers = ""
        def characters = ""
        last_digit.each { d ->
            if(d.isNumber()) {
                numbers += d
            } else {
                characters += d
            }
        }
        init_list[-1] = numbers
        init_list.add(characters)
    }
    return init_list.collect { num -> num.isNumber() ? num.toInteger() : num }
}
