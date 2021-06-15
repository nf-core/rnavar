/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.


def valid_params = [
    aligners       : ['star'],
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnavar.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, 
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.star_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner }

def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]


def star_genomegenerate_options = modules['star_genomegenerate']
if (!params.save_reference)     { star_genomegenerate_options['publish_files'] = false }

def samtools_sort_genome_options  = modules['samtools_sort_genome']
def samtools_index_genome_options = modules['samtools_index_genome']
samtools_index_genome_options.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''
if (['star'].contains(params.aligner)) {
    if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
        samtools_sort_genome_options.publish_files.put('bam','')
        samtools_index_genome_options.publish_files.put('bai','')
        samtools_index_genome_options.publish_files.put('csi','')
    }
} else {
    if (params.save_align_intermeds || params.skip_markduplicates) {
        samtools_sort_genome_options.publish_files.put('bam','')
        samtools_index_genome_options.publish_files.put('bai','')
        samtools_index_genome_options.publish_files.put('csi','')
    }
}

def biotype                       = params.gencode ? "gene_type" : params.featurecounts_group_type

def cat_fastq_options             = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

def star_align_options            = modules['star_align']
star_align_options.args          += params.save_unaligned ? Utils.joinModuleArgs(['--outReadsUnmapped Fastx']) : ''
if (params.save_align_intermeds)  { star_align_options.publish_files.put('bam','') }
if (params.save_unaligned)        { star_align_options.publish_files.put('fastq.gz','unmapped') }

def picard_markduplicates_samtools   = modules['picard_markduplicates_samtools']
picard_markduplicates_samtools.args += params.bam_csi_index ? Utils.joinModuleArgs(['-c']) : ''

include { INPUT_CHECK } from '../subworkflows/local/input_check'    addParams( options: [:] )
include { CAT_FASTQ             } from '../modules/nf-core/software/cat/fastq/main' addParams( options: cat_fastq_options )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams( genome_options: publish_genome_options, index_options: publish_index_options, star_index_options: star_genomegenerate_options)
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/software/gatk4/bedtointervallist/main' addParams( options: [:])
include { SAMTOOLS_INDEX } from '../modules/nf-core/software/samtools/index/main' addParams( options: samtools_index_genome_options )
include { ALIGN_STAR } from '../subworkflows/nf-core/align_star'       addParams( align_options: star_align_options, samtools_sort_options: samtools_sort_genome_options, samtools_index_options: samtools_index_genome_options, samtools_stats_options: samtools_index_genome_options   )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: picard_markduplicates_samtools, samtools_stats_options:  picard_markduplicates_samtools )
include { GATK4_SPLITNCIGARREADS } from '../modules/nf-core/software/gatk4/splitncigarreads/main' addParams( options: publish_genome_options, splitncigar_options: modules['gatk_splitncigar_options'] )
include { GATK4_BASERECALIBRATOR } from '../modules/nf-core/software/gatk4/baserecalibrator/main' addParams( options: modules['gatk_baserecalibrator'])

workflow RNASEQ_VAR {

    PREPARE_GENOME (
        prepareToolIndices,
        biotype
    )
    ch_genome_gtf = Channel.empty()
    ch_genome_gtf = PREPARE_GENOME.out.gene_bed
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))


    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq } 

    // PREPARE THE INTERVAL LIST
    ch_interval_list = Channel.empty()
    GATK4_BEDTOINTERVALLIST(
        ch_genome_gtf,
        PREPARE_GENOME.out.dict
    )
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_software_versions      = ch_software_versions.mix(GATK4_BEDTOINTERVALLIST.out.version.first().ifEmpty(null))
    //
    // SUBWORKFLOW: Alignment with STAR
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()

    if (!params.skip_alignment && params.aligner == 'star') {
        ALIGN_STAR (
            ch_cat_fastq,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        if (params.bam_csi_index) {
            ch_genome_bam_index  = ALIGN_STAR.out.csi
        }
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))

        // SUBWORKFLOW: Mark duplicates with Picard
        ch_markduplicates_multiqc = Channel.empty()
        if (!params.skip_alignment && !params.skip_markduplicates) {
            MARK_DUPLICATES_PICARD (
                ch_genome_bam
            )
            ch_genome_bam             = MARK_DUPLICATES_PICARD.out.bam
            ch_genome_bam_index       = MARK_DUPLICATES_PICARD.out.bai
            ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
            ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
            ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
            ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
            if (params.bam_csi_index) {
                ch_genome_bam_index  = MARK_DUPLICATES_PICARD.out.csi
            }
            ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
        }
        // MODULE: SplitNCigarReads from GATK4
        // Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data). 

        if(!params.skip_splitncigar){

            GATK4_SPLITNCIGARREADS (
                ch_genome_bam,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict
            )

            SAMTOOLS_INDEX ( GATK4_SPLITNCIGARREADS.out.bam )

            ch_genome_bam           = GATK4_SPLITNCIGARREADS.out.bam
            ch_genome_bam_index     = SAMTOOLS_INDEX.out.bai
            ch_software_versions    = ch_software_versions.mix(GATK4_SPLITNCIGARREADS.out.version.first().ifEmpty(null))
        }

        /*
        if(!params.skip_basecalibrator){
            GATK4_BASERECALIBRATOR(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fai,
                PREPARE_GENOME.out.dict,
                ch_interval_list,
            )
        }*/
    }
}