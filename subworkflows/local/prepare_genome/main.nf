//
// Prepare reference genome files
//

include { GATK4_CREATESEQUENCEDICTIONARY                       } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GFFREAD                                              } from '../../../modules/nf-core/gffread'
include { GTF2BED                                              } from '../../../modules/local/gtf2bed'
include { GUNZIP as GUNZIP_FASTA                               } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                                 } from '../../../modules/nf-core/gunzip'
include { HTSLIB_BGZIPTABIX as BGZIPTABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/htslib/bgziptabix'
include { HTSLIB_BGZIPTABIX as BGZIPTABIX_DBSNP                } from '../../../modules/nf-core/htslib/bgziptabix'
include { HTSLIB_BGZIPTABIX as BGZIPTABIX_KNOWN_INDELS         } from '../../../modules/nf-core/htslib/bgziptabix'
include { REMOVEUNKNOWNREGIONS                                 } from '../../../modules/local/removeunknownregions'
include { SAMTOOLS_FAIDX                                       } from '../../../modules/nf-core/samtools/faidx'
include { STAR_GENOMEGENERATE                                  } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_INDEXVERSION                                    } from '../../../modules/nf-core/star/indexversion'
include { UNTAR                                                } from '../../../modules/nf-core/untar'

workflow PREPARE_GENOME {
    take:
    bcftools_annotations // params[path]: params.bcftools_annotations
    bcftools_annotations_tbi // params[path]: params.bcftools_annotations_tbi
    dbsnp // params[path]: params.dbsnp
    dbsnp_tbi // params[path]: params.dbsnp_tbi
    dict // params[path]: params.dict
    exon_bed // params[path]: params.exon_bed
    fasta // params[path]: params.fasta
    fasta_fai // params[path]: params.fasta_fai
    gff // params[path]: params.gff
    gtf // params[path]: params.gtf
    known_indels // params[path]: params.known_indels
    known_indels_tbi // params[path]: params.known_indels_tbi
    star_index // params[path]: params.star_index
    feature_type // params[val]: params.feature_type
    align // boolean: The pipeline needs aligner indices or not
    genome // params[val]: params.genome
    tools // list: list of tools to run

    main:
    // Unzip reference genome files if needed
    def ch_gunzip_fasta_input = fasta.toString().endsWith('.gz')
        ? channel.fromPath(fasta).map { fasta_ -> [[id: genome], fasta_] }.collect()
        : channel.empty()

    GUNZIP_FASTA(ch_gunzip_fasta_input)

    def ch_fasta = fasta.toString().endsWith('.gz')
        ? GUNZIP_FASTA.out.gunzip.collect()
        : channel.fromPath(fasta).map { fasta_ -> [[id: genome], fasta_] }.collect()

    def dict_input = dict ? channel.empty() : ch_fasta

    GATK4_CREATESEQUENCEDICTIONARY(dict_input)

    def ch_dict = dict
        ? channel.fromPath(dict).map { dict_ -> [[id: genome], dict_] }.collect()
        : GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()

    def gtf_input = gtf.toString().endsWith('.gz')
        ? channel.fromPath(gtf).map { gtf_ -> [[id: genome], gtf_] }
        : channel.empty()

    GUNZIP_GTF(gtf_input)

    def ch_gffread_input = gff
        ? channel.fromPath(gff).map { gff_ -> [[id: genome], gff_] }
        : channel.empty()

    GFFREAD(ch_gffread_input, ch_fasta.map { _meta, fasta_ -> fasta_ })

    def ch_gtf = gtf.toString().endsWith('.gz')
        ? GUNZIP_GTF.out.gunzip.collect()
        : gff
            ? GFFREAD.out.gtf.collect()
            : channel.fromPath(gtf).map { gtf_ -> [[id: genome], gtf_] }.collect()

    def ch_gtf2bed_input = !exon_bed ? ch_gtf : channel.empty()

    GTF2BED(ch_gtf2bed_input, feature_type)

    def ch_exon_bed_input = exon_bed
        ? channel.fromPath(exon_bed).map { exon_bed_ -> [[id: genome], exon_bed_] }.collect()
        : GTF2BED.out.bed.collect()

    def ch_remove_unknown_regions_input = 'removeunknownregions' in tools ? ch_exon_bed_input : channel.empty()

    REMOVEUNKNOWNREGIONS(ch_remove_unknown_regions_input.join(ch_dict))

    def ch_exon_bed = 'removeunknownregions' in tools ? REMOVEUNKNOWNREGIONS.out.bed : ch_exon_bed_input

    def ch_bcftools_annotations_in = bcftools_annotations
        ? channel.fromPath(bcftools_annotations)
        : channel.value([])
    def ch_bcftools_annotations_tbi = bcftools_annotations_tbi
        ? channel.fromPath(bcftools_annotations_tbi).collect()
        : channel.value([])

    // we use vcf.baseName - '.vcf', because we have to deal with both .vcf and .vcf.gz
    if (!bcftools_annotations_tbi && bcftools_annotations) {
        BGZIPTABIX_BCFTOOLS_ANNOTATIONS(
            ch_bcftools_annotations_in.map { vcf -> [[id: vcf.baseName - '.vcf'], vcf, [], []] },
            'compress',
            true,
            'vcf',
        )
        ch_bcftools_annotations_tbi = BGZIPTABIX_BCFTOOLS_ANNOTATIONS.out.index.map { _meta, tbi -> [tbi] }.collect()
        ch_bcftools_annotations_vcf = BGZIPTABIX_BCFTOOLS_ANNOTATIONS.out.output.map { _meta, vcf -> [vcf] }.collect()
    }
    else {
        ch_bcftools_annotations_vcf = ch_bcftools_annotations_in.collect()
    }

    // we use vcf.baseName - '.vcf', because we have to deal with both .vcf and .vcf.gz
    def ch_dbsnp_in = dbsnp
        ? channel.fromPath(dbsnp).flatten().map { vcf -> [[id: vcf.baseName - '.vcf'], vcf] }
        : channel.value([[id: genome], []])
    def ch_dbsnp_tbi = dbsnp_tbi
        ? channel.fromPath(dbsnp_tbi).flatten().map { tbi -> [[id: genome], tbi] }
        : channel.value([[id: genome], []])

    if (!dbsnp_tbi && dbsnp) {
        BGZIPTABIX_DBSNP(
            ch_dbsnp_in.map { meta, vcf -> [meta, vcf, [], []] },
            'compress',
            true,
            'vcf',
        )
        ch_dbsnp_tbi = BGZIPTABIX_DBSNP.out.index
        ch_dbsnp_vcf = BGZIPTABIX_DBSNP.out.output.map { meta, file -> [meta + [id: genome], file] }
    }
    else {
        ch_dbsnp_vcf = ch_dbsnp_in.map { meta, file -> [meta + [id: genome], file] }
    }

    // we use vcf.baseName - '.vcf', because we have to deal with both .vcf and .vcf.gz
    def ch_known_indels_in = known_indels
        ? channel.fromPath(known_indels).flatten().map { vcf -> [[id: vcf.baseName - '.vcf'], vcf] }
        : channel.value([[id: genome], []])
    def ch_known_indels_tbi = known_indels_tbi
        ? channel.fromPath(known_indels_tbi).flatten().map { tbi -> [[genome], tbi] }
        : channel.value([[id: genome], []])

    if (!known_indels_tbi && known_indels) {
        BGZIPTABIX_KNOWN_INDELS(
            ch_known_indels_in.map { meta, vcf -> [meta, vcf, [], []] },
            'compress',
            true,
            'vcf',
        )
        ch_known_indels_tbi = BGZIPTABIX_KNOWN_INDELS.out.index
        ch_known_indels_vcf = BGZIPTABIX_KNOWN_INDELS.out.output.map { meta, file -> [meta + [id: genome], file] }
    }
    else {
        ch_known_indels_vcf = ch_known_indels_in.map { meta, file -> [meta + [id: genome], file] }
    }

    // known_sites is made by grouping both the dbsnp and the known indels resources
    // Which can either or both be optional
    def ch_known_sites_vcf = ch_dbsnp_vcf
        .mix(ch_known_indels_vcf)
        .collect { _meta, file -> file }
        .map { file -> [[id: genome], file] }

    def ch_known_sites_tbi = ch_dbsnp_tbi
        .mix(ch_known_indels_tbi)
        .collect { _meta, file -> file }
        .map { file -> [[id: genome], file] }

    def fai_input = fasta_fai
        ? channel.empty()
        : ch_fasta.map { meta, _fasta -> [meta, _fasta, []] }

    SAMTOOLS_FAIDX(fai_input, false)

    def ch_fai = fasta_fai
        ? channel.fromPath(fasta_fai).map { fai_ -> [[id: genome], fai_] }.collect()
        : SAMTOOLS_FAIDX.out.fai.collect()

    //
    // STAR index handling
    //

    def star_index_input = star_index
        ? channel.fromPath(star_index).map { index -> [[id: genome], index] }
        : channel.of([[], []])

    ch_star_index_input = star_index_input
        .map { _meta, index -> [[id: genome], index] }
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

    UNTAR(ch_star_index_input.tarzipped)

    STAR_INDEXVERSION()

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

    star_index_output = STAR_GENOMEGENERATE.out.index
        .mix(star_index_check.compatible)
        .collect()

    emit:
    bcfann           = ch_bcftools_annotations_vcf // path: bcftools_annotations.vcf.gz
    bcfann_tbi       = ch_bcftools_annotations_tbi // path: bcftools_annotations.vcf.gz.tbi
    dbsnp            = ch_dbsnp_vcf.collect() // Channel: [meta, dbsnp.vcf.gz]
    dbsnp_tbi        = ch_dbsnp_tbi.collect() // Channel: [meta, dbsnp.vcf.gz.tbi]
    dict             = ch_dict // path: genome.fasta.dict
    exon_bed         = ch_exon_bed // path: exon.bed
    fasta            = ch_fasta // path: genome.fasta
    fasta_fai        = ch_fai // path: genome.fasta.fai
    gtf              = ch_gtf // path: genome.gtf
    known_indels     = ch_known_indels_vcf.collect() // path: {known_indels*}.vcf.gz
    known_indels_tbi = ch_known_indels_tbi.collect() // path: {known_indels*}.vcf.gz.tbi
    known_sites      = ch_known_sites_vcf // path: {known_sites*}.vcf.gz
    known_sites_tbi  = ch_known_sites_tbi // path: {known_sites*}.vcf.gz.tbi
    star_index       = star_index_output // path: star/index/
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
