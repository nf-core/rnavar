//
// Prepare reference genome files
//

include { GATK4_CREATESEQUENCEDICTIONARY                      } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GFFREAD                                             } from '../../../modules/nf-core/gffread'
include { GTF2BED                                             } from '../../../modules/local/gtf2bed'
include { GUNZIP as GUNZIP_FASTA                              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                                } from '../../../modules/nf-core/gunzip'
include { REMOVE_UNKNOWN_REGIONS                              } from '../../../modules/local/remove_unkown_regions'
include { SAMTOOLS_FAIDX                                      } from '../../../modules/nf-core/samtools/faidx'
include { STAR_GENOMEGENERATE                                 } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_INDEXVERSION                                   } from '../../../modules/nf-core/star/indexversion'
include { TABIX_BGZIPTABIX as BGZIPTABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_DBSNP                } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_KNOWN_INDELS         } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX as TABIX_BCFTOOLS_ANNOTATIONS           } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_DBSNP                          } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_KNOWN_INDELS                   } from '../../../modules/nf-core/tabix/tabix'
include { UNTAR                                               } from '../../../modules/nf-core/untar'

workflow PREPARE_GENOME {
    take:
    fasta                    // params[path]: params.fasta
    dict                     // params[path]: params.dict
    fai                      // params[path]: params.fasta_fai
    star_index               // params[path]: params.star_index
    gff                      // params[path]: params.gff
    gtf                      // params[path]: params.gtf
    exon_bed                 // params[path]: params.exon_bed
    bcftools_annotations     // params[path]: params.bcftools_annotations
    bcftools_annotations_tbi // params[path]: params.bcftools_annotations_tbi
    dbsnp                    // params[path]: params.dbsnp
    dbsnp_tbi                // params[path]: params.dbsnp_tbi
    known_indels             // params[path]: params.known_indels
    known_indels_tbi         // params[path]: params.known_indels_tbi
    feature_type             // params[string]: params.feature_type
    skip_exon_bed_check      // params[boolean]: params.skip_exon_bed_check
    align                    // boolean: The pipeline needs aligner indices or not

    main:
    def ch_versions = Channel.empty()

    //Unzip reference genome files if needed

    def ch_fasta = Channel.empty()
    if (fasta.toString().endsWith('.gz')) {
        GUNZIP_FASTA(
            fasta.map { fasta_ -> [[id: fasta_.baseName], fasta_] }
        )

        ch_fasta = GUNZIP_FASTA.out.gunzip.collect()
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }
    else {
        ch_fasta = fasta
            .map { fasta_ -> [[id: fasta_.baseName], fasta_] }
            .collect()
    }

    def ch_dict = Channel.empty()
    if (!dict) {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)

        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }
    else {
        ch_dict = dict
            .map { dict_ -> [[id: dict_.baseName], dict_] }
            .collect()
    }

    def ch_gtf = Channel.empty()
    if (gtf.toString().endsWith('.gz')) {
        GUNZIP_GTF(
            gtf.map { gtf_ -> [[id: gtf_.baseName], gtf_] }
        )

        ch_gtf = GUNZIP_GTF.out.gunzip.collect()
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    }
    else if (gff) {
        GFFREAD(
            gff.map { gff_ -> [[id: gff_.baseName], gff_] },
            ch_fasta.map { _meta, fasta_ -> fasta_ },
        )

        ch_gtf = GFFREAD.out.gtf.collect()
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }
    else {
        ch_gtf = gtf
            .map { gtf_ -> [[id: gtf_.baseName], gtf_] }
            .collect()
    }

    def ch_exon_bed_raw = Channel.empty()
    if (!exon_bed) {
        GTF2BED(ch_gtf, feature_type)

        ch_exon_bed_raw = GTF2BED.out.bed.collect()
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }
    else {
        ch_exon_bed_raw = exon_bed
            .map { exon_bed_ -> [[id: exon_bed_.baseName], exon_bed_] }
            .collect()
    }

    def ch_exon_bed = Channel.empty()
    if (!skip_exon_bed_check) {
        REMOVE_UNKNOWN_REGIONS(
            ch_exon_bed_raw,
            ch_dict,
        )

        ch_exon_bed = REMOVE_UNKNOWN_REGIONS.out.bed.collect()
        ch_versions = ch_versions.mix(REMOVE_UNKNOWN_REGIONS.out.versions)
    }
    else {
        ch_exon_bed = ch_exon_bed_raw
    }

    def bcftools_annotations_input = bcftools_annotations
        ? bcftools_annotations.map { vcf -> [[id: vcf.name], vcf] }
        : Channel.empty()

    def bcftools_annotations_tbi_input = bcftools_annotations_tbi
        ? bcftools_annotations_tbi.map { tbi -> [[id: tbi.baseName], tbi] }
        : Channel.empty()

    def ch_bcftools_annotations_input = bcftools_annotations_input
        .join(bcftools_annotations_tbi_input, failOnDuplicate: true, remainder: true)
        .branch { meta, file, index ->
            plain: !file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_noindex: !index && file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_index: true
            return [meta, file, index]
        }

    BGZIPTABIX_BCFTOOLS_ANNOTATIONS(
        ch_bcftools_annotations_input.plain
    )
    ch_versions = ch_versions.mix(BGZIPTABIX_BCFTOOLS_ANNOTATIONS.out.versions)

    TABIX_BCFTOOLS_ANNOTATIONS(
        ch_bcftools_annotations_input.bgzip_noindex
    )
    ch_versions = ch_versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)

    def ch_bcftools_annotations = BGZIPTABIX_BCFTOOLS_ANNOTATIONS.out.gz_tbi
        .map { _meta, file, _index -> file }
        .mix(ch_bcftools_annotations_input.bgzip_noindex.map { _meta, file -> file })
        .mix(ch_bcftools_annotations_input.bgzip_index.map { _meta, file, _index -> file })
        .collect()

    def ch_bcftools_annotations_tbi = BGZIPTABIX_BCFTOOLS_ANNOTATIONS.out.gz_tbi
        .map { _meta, _file, index -> index }
        .mix(TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map { _meta, tbi -> tbi })
        .mix(ch_bcftools_annotations_input.bgzip_index.map { _meta, _file, index -> index })
        .collect()

    def dbsnp_input = dbsnp
        ? dbsnp.flatten().map { vcf -> [[id: vcf.name], vcf] }
        : Channel.empty()

    def dbsnp_tbi_input = dbsnp_tbi
        ? dbsnp_tbi.flatten().map { tbi -> [[id: tbi.baseName], tbi] }
        : Channel.empty()

    def ch_dbsnp_input = dbsnp_input
        .join(dbsnp_tbi_input, failOnDuplicate: true, remainder: true)
        .branch { meta, file, index ->
            plain: !file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_noindex: !index && file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_index: true
            return [meta, file, index]
        }

    BGZIPTABIX_DBSNP(
        ch_dbsnp_input.plain
    )
    ch_versions = ch_versions.mix(BGZIPTABIX_DBSNP.out.versions)

    TABIX_DBSNP(
        ch_dbsnp_input.bgzip_noindex
    )
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)

    def ch_dbsnp = BGZIPTABIX_DBSNP.out.gz_tbi
        .map { _meta, file, _index -> file }
        .mix(ch_dbsnp_input.bgzip_noindex.map { _meta, file -> file })
        .mix(ch_dbsnp_input.bgzip_index.map { _meta, file, _index -> file })
        .collect()

    def ch_dbsnp_tbi = BGZIPTABIX_DBSNP.out.gz_tbi
        .map { _meta, _file, index -> index }
        .mix(TABIX_DBSNP.out.tbi.map { _meta, tbi -> tbi })
        .mix(ch_dbsnp_input.bgzip_index.map { _meta, _file, index -> index })
        .collect()

    def known_indels_input = known_indels
        ? known_indels.flatten().map { vcf -> [[id: vcf.name], vcf] }
        : Channel.empty()

    def known_indels_tbi_input = known_indels_tbi
        ? known_indels_tbi.flatten().map { tbi -> [[id: tbi.baseName], tbi] }
        : Channel.empty()

    def ch_known_indels_input = known_indels_input
        .join(known_indels_tbi_input, failOnDuplicate: true, remainder: true)
        .branch { meta, file, index ->
            plain: !file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_noindex: !index && file.toString().endsWith(".gz")
            return [meta, file]
            bgzip_index: true
            return [meta, file, index]
        }

    BGZIPTABIX_KNOWN_INDELS(
        ch_known_indels_input.plain
    )
    ch_versions = ch_versions.mix(BGZIPTABIX_KNOWN_INDELS.out.versions)

    TABIX_KNOWN_INDELS(
        ch_known_indels_input.bgzip_noindex
    )
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)

    def ch_known_indels = BGZIPTABIX_KNOWN_INDELS.out.gz_tbi
        .map { _meta, file, _index -> file }
        .mix(ch_known_indels_input.bgzip_noindex.map { _meta, file -> file })
        .mix(ch_known_indels_input.bgzip_index.map { _meta, file, _index -> file })
        .collect()

    def ch_known_indels_tbi = BGZIPTABIX_KNOWN_INDELS.out.gz_tbi
        .map { _meta, _file, index -> index }
        .mix(TABIX_KNOWN_INDELS.out.tbi.map { _meta, tbi -> tbi })
        .mix(ch_known_indels_input.bgzip_index.map { _meta, _file, index -> index })
        .collect()

    def ch_fai = Channel.empty()
    if (!fai) {
        SAMTOOLS_FAIDX(ch_fasta, [[id: ch_fasta.baseName], []], false)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai = SAMTOOLS_FAIDX.out.fai.collect()
    }
    else {
        ch_fai = fai
            .map { fai_ -> [[id: fai_.baseName], fai_] }
            .collect()
    }

    //
    // STAR index handling
    //

    def star_index_input = star_index
        ? star_index.map { index -> [[id: 'star'], index] }
        : Channel.of([[], []])

    ch_star_index_input = star_index_input
        .map { _meta, index -> [[id: 'star'], index] }
        .merge(align)
        .filter { _meta, _index, bool_align ->
            return bool_align
        }
        .branch { meta, index, _align ->
            no_index: !index
            return [meta, index]
            tarzipped: index.name.endsWith(".tar.gz")
            return [meta, index]
            index: true
            return [meta, index]
        }

    UNTAR(
        ch_star_index_input.tarzipped
    )
    ch_versions = ch_versions.mix(UNTAR.out.versions)

    STAR_INDEXVERSION()
    ch_versions = ch_versions.mix(STAR_INDEXVERSION.out.versions)

    def star_index_check = ch_star_index_input.index
        .mix(UNTAR.out.untar)
        .combine(STAR_INDEXVERSION.out.index_version)
        .branch { meta, index, version_file ->
            def is_compatible = true
            if (!workflow.stubRun) {
                def minimal_version = version_file.text.replace("\n", "")
                def index_version = index.resolve("genomeParameters.txt").text.readLines().find { line -> line.startsWith("versionGenome") }.tokenize("\t")[-1]
                is_compatible = isCompatibleStarIndex(index_version, minimal_version)
                if (!is_compatible) {
                    log.warn("Detected a wrong version of the STAR index, expected a minimum version of ${minimal_version}. Automatically recreating the index of STAR...")
                }
            }
            compatible: is_compatible
            return [meta, index]
            incompatible: !is_compatible
            return [meta, []]
        }

    def genomegenerate_input = star_index_check.incompatible
        .mix(ch_star_index_input.no_index)
        .combine(ch_fasta)
        .map { _meta1, _wrong_index, meta2, fasta_ ->
            [meta2, fasta_]
        }

    STAR_GENOMEGENERATE(genomegenerate_input, ch_gtf)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    star_index_output = STAR_GENOMEGENERATE.out.index
        .mix(star_index_check.compatible)
        .collect()

    emit:
    bcfann           = ch_bcftools_annotations // path: bcftools_annotations.vcf.gz
    bcfann_tbi       = ch_bcftools_annotations_tbi // path: bcftools_annotations.vcf.gz.tbi
    dbsnp            = ch_dbsnp // path: dbsnp.vcf.gz
    dbsnp_tbi        = ch_dbsnp_tbi // path: dbsnp.vcf.gz.tbi
    dict             = ch_dict // path: genome.fasta.dict
    exon_bed         = ch_exon_bed // path: exon.bed
    fasta            = ch_fasta // path: genome.fasta
    fasta_fai        = ch_fai // path: genome.fasta.fai
    gtf              = ch_gtf // path: genome.gtf
    known_indels     = ch_known_indels // path: {known_indels*}.vcf.gz
    known_indels_tbi = ch_known_indels_tbi // path: {known_indels*}.vcf.gz.tbi
    star_index       = star_index_output // path: star/index/
    versions         = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if the STAR index is compatible with the minimal version
def isCompatibleStarIndex(index_version, minimal_index_version) {
    def is_compatible = true
    if (minimal_index_version.isNumber()) {
        // Older version of STAR used a numerical versioning.
        // Return true if the index doesn't use the numerical versioning anymore
        if (!index_version.isNumber()) {
            is_compatible = true
        }
        else {
            is_compatible = index_version.toInteger() >= minimal_index_version.toInteger()
        }
    }
    else {
        if (index_version.isNumber()) {
            is_compatible = false
        }
        else {
            // Correctly compare semantic version strings: e.g 2.7.11b > 2.7.4a
            def min_list = convertVersionToList(minimal_index_version)
            def ind_list = convertVersionToList(index_version)
            ind_list.eachWithIndex { digit, idx ->
                if (digit > min_list[idx]) {
                    is_compatible = true
                    return null
                }
                else if (digit < min_list[idx]) {
                    is_compatible = false
                    return null
                }
            }
        }
    }
    return is_compatible
}

// Convert a version string to a list of numbers and characters
def convertVersionToList(version) {
    def init_list = version.tokenize(".")
    if (!init_list[-1].isNumber()) {
        // Handle cases where the last digit in the version contains a character: e.g. 2.7.11b
        def last_digit = init_list[-1]
        def numbers = ""
        def characters = ""
        last_digit.each { d ->
            if (d.isNumber()) {
                numbers += d
            }
            else {
                characters += d
            }
        }
        init_list[-1] = numbers
        init_list.add(characters)
    }
    return init_list.collect { num -> num.isNumber() ? num.toInteger() : num }
}
