/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { GTF2BED                   } from '../../modules/local/gtf2bed'

// nf-core
include { CAT_FASTQ                 } from '../../modules/nf-core/cat/fastq'
include { FASTQC                    } from '../../modules/nf-core/fastqc'
include { GATK4_BASERECALIBRATOR    } from '../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_BEDTOINTERVALLIST   } from '../../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_COMBINEGVCFS        } from '../../modules/nf-core/gatk4/combinegvcfs'
include { GATK4_HAPLOTYPECALLER     } from '../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_INDEXFEATUREFILE    } from '../../modules/nf-core/gatk4/indexfeaturefile'
include { GATK4_INTERVALLISTTOOLS   } from '../../modules/nf-core/gatk4/intervallisttools'
include { GATK4_MERGEVCFS           } from '../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_VARIANTFILTRATION   } from '../../modules/nf-core/gatk4/variantfiltration'
include { MULTIQC                   } from '../../modules/nf-core/multiqc'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index'
include { TABIX_TABIX as TABIX      } from '../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIXGVCF  } from '../../modules/nf-core/tabix/tabix'

// local
include { RECALIBRATE               } from '../../subworkflows/local/recalibrate'
include { SPLITNCIGAR               } from '../../subworkflows/local/splitncigar'
include { VCF_ANNOTATE_ALL          } from '../../subworkflows/local/vcf_annotate_all'
include { PREPARE_ALIGNMENT         } from '../../subworkflows/local/prepare_alignment'

// nf-core
include { BAM_MARKDUPLICATES_PICARD } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { FASTQ_ALIGN_STAR          } from '../../subworkflows/nf-core/fastq_align_star'

// local
include { checkSamplesAfterGrouping } from '../../subworkflows/local/utils_nfcore_rnavar_pipeline'
include { methodsDescriptionText    } from '../../subworkflows/local/utils_nfcore_rnavar_pipeline'

// nf-core
include { paramsSummaryMultiqc      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

// plugin
include { paramsSummaryMap          } from 'plugin/nf-schema'

/*
========================================================================================
    RUN MAIN WORKFLOW RNAVAR
========================================================================================
*/

workflow RNAVAR {
    take:
    input
    dbsnp
    dbsnp_tbi
    dict
    exon_bed
    fasta
    fasta_fai
    gtf
    known_indels
    known_indels_tbi
    star_index
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files
    seq_center
    seq_platform

    main:

    // To gather all QC reports for MultiQC
    def ch_reports  = Channel.empty()

    // To gather used softwares ch_versions for MultiQC
    def ch_versions = Channel.empty()

    // Parse the input data
    parsed_input = input.groupTuple().map{ samplesheet -> checkSamplesAfterGrouping(samplesheet) }
        .branch{ meta, fastqs, bam, bai, cram, crai ->
            single  : fastqs.size() == 1
                return [ meta, fastqs.flatten() ]
            multiple: fastqs.size() > 1
                return [ meta, fastqs.flatten() ]
            bam     : bam
                return [ meta, bam, bai ]
            cram    : cram
                return [ meta, cram, crai ]
        }

    // MODULE: Prepare the alignment files (convert CRAM -> BAM and index)
    PREPARE_ALIGNMENT(
        parsed_input.cram,
        parsed_input.bam,
    )
    ch_versions = ch_versions.mix(PREPARE_ALIGNMENT.out.versions)

    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ(parsed_input.multiple)

    def cat_fastq = CAT_FASTQ.out.reads.mix(parsed_input.single)

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // MODULE: Generate QC summary using FastQC
    FASTQC(cat_fastq)
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{ _meta, logs -> logs })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
    //

    GATK4_BEDTOINTERVALLIST(exon_bed, dict)
    def interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    //
    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    //
    def interval_list_split = Channel.empty()
    if (!params.skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(interval_list)
        interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ _meta, bed -> [bed] }.collect()
    }
    else {
        interval_list_split = interval_list.map { _meta, bed -> bed }
    }

    //
    // SUBWORKFLOW: Perform read alignment using STAR aligner
    //

    if (params.aligner == 'star') {
        FASTQ_ALIGN_STAR(cat_fastq,
            star_index,
            gtf,
            params.star_ignore_sjdbgtf,
            seq_platform,
            seq_center,
            fasta,
            [[:],[]]) //transcripts_fasta)

        def genome_bam    = FASTQ_ALIGN_STAR.out.bam

        // Gather QC ch_reports
        ch_reports           = ch_reports.mix(FASTQ_ALIGN_STAR.out.log_out)
        ch_reports           = ch_reports.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_versions          = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Mark duplicates with GATK4
        //
        BAM_MARKDUPLICATES_PICARD(genome_bam,
            fasta,
            fasta_fai)

        def markduplicate_indices = BAM_MARKDUPLICATES_PICARD.out.bai
            .mix(BAM_MARKDUPLICATES_PICARD.out.csi)
            .mix(BAM_MARKDUPLICATES_PICARD.out.crai)

        def genome_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(markduplicate_indices, failOnDuplicate:true, failOnMismatch:true)
            .mix(PREPARE_ALIGNMENT.out.bam)

        //Gather QC ch_reports
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]))
        ch_versions               = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        //
        // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string(e.g. spanning splicing events in RNAseq data).
        //

        SPLITNCIGAR(genome_bam_bai,
            fasta,
            fasta_fai.map { _meta, fai -> fai },
            dict,
            interval_list.map { _meta, bed -> bed }
        )

        def splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
        ch_versions                 = ch_versions.mix(SPLITNCIGAR.out.versions)

        //
        // MODULE: BaseRecalibrator from GATK4
        // Generates a recalibration table based on various co-variates
        //
        def bam_variant_calling = Channel.empty()

        if (!params.skip_baserecalibration) {
            // known_sites is made by grouping both the dbsnp and the known indels ressources
            // they can either or both be optional
            def known_sites     = dbsnp.concat(known_indels).collect()
            def known_sites_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

            def interval_list_recalib = interval_list.map{ _meta, bed -> [bed] }.flatten()
            def splitncigar_bam_bai_interval = splitncigar_bam_bai.combine(interval_list_recalib)
                .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval] }

            GATK4_BASERECALIBRATOR(
                splitncigar_bam_bai_interval,
                fasta.map{ _meta, fasta_ -> [fasta_] },
                fasta_fai.map { _meta, fai -> fai },
                dict.map{ _meta, dict_ -> [dict_] },
                known_sites,
                known_sites_tbi
            )
            def bqsr_table   = GATK4_BASERECALIBRATOR.out.table

            // Gather QC ch_reports
            ch_reports  = ch_reports.mix(bqsr_table.map{ _meta, table -> table})
            ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

            def bam_applybqsr       = splitncigar_bam_bai.join(bqsr_table)

            def interval_list_applybqsr = interval_list.map{ _meta, bed -> [bed] }.flatten()
            def applybqsr_bam_bai_interval = bam_applybqsr.combine(interval_list_applybqsr)
                .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]}

            //
            // MODULE: ApplyBaseRecalibrator from GATK4
            // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
            //
            RECALIBRATE(
                params.skip_multiqc,
                applybqsr_bam_bai_interval,
                dict.map{ _meta, dict_ -> [dict_] },
                fasta_fai.map { _meta, fai -> fai },
                fasta.map{ _meta, fasta_ -> [fasta_] }
            )

            bam_variant_calling = RECALIBRATE.out.bam

            // Gather QC ch_reports
            ch_reports  = ch_reports.mix(RECALIBRATE.out.qc.collect{it[1]}.ifEmpty([]))
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions)
        } else {
            bam_variant_calling = splitncigar_bam_bai
        }

        // Run haplotyper even in the absence of dbSNP files
        def dbsnp_for_haplotypecaller = Channel.empty()
        def dbsnp_for_haplotypecaller_tbi = Channel.empty()
        if (!params.dbsnp){
            dbsnp_for_haplotypecaller = [[id:'null'], []]
            dbsnp_for_haplotypecaller_tbi = [[id:'null'], []]
        } else {
            dbsnp_for_haplotypecaller     = dbsnp.map{ vcf -> [[id:'dbsnp'], vcf] }
            dbsnp_for_haplotypecaller_tbi = dbsnp_tbi.map{ tbi -> [[id:'dbsnp'], tbi] }
        }

        def haplotypecaller_interval_bam = bam_variant_calling.combine(interval_list_split)
            .map { meta, bam, bai, interval_lists ->
                def new_meta = meta + [interval_count: interval_lists instanceof List ? interval_lists.size() : 1]
                [ new_meta, bam, bai, interval_lists ]
            }
            .transpose(by:3)
            .map{ meta, bam, bai, interval_list_ ->
                [ meta + [ id:meta.id + "_" + interval_list_.baseName, sample:meta.id, variantcaller:'haplotypecaller' ], bam, bai, interval_list_, [] ]
            }

        //
        // MODULE: HaplotypeCaller from GATK4
        // Calls germline SNPs and indels via local re-assembly of haplotypes.
        //

        GATK4_HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            fasta,
            fasta_fai,
            dict,
            dbsnp_for_haplotypecaller,
            dbsnp_for_haplotypecaller_tbi
        )

        def haplotypecaller_out = GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .map{ meta, vcf, tbi -> 
                def new_meta = meta + [id:meta.sample] - meta.subMap('sample') - meta.subMap("interval_count")
                [ groupKey(new_meta, meta.interval_count), vcf, tbi ]
            }
            .groupTuple()

        ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

        def haplotypecaller_vcf = Channel.empty()
        if (!params.generate_gvcf){
            //
            // MODULE: MergeVCFS from GATK4
            // Merge multiple VCF files into one VCF
            //
            def haplotypecaller_raw = haplotypecaller_out.map { meta, vcfs, _tbis -> [ meta, vcfs ]}
            GATK4_MERGEVCFS(
                haplotypecaller_raw,
                dict
            )
            haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
            ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

            //
            // MODULE: Index the VCF using TABIX
            //
            TABIX(
                haplotypecaller_vcf
            )
            ch_versions      = ch_versions.mix(TABIX.out.versions)

            def haplotypecaller_indices = TABIX.out.tbi
                .mix(TABIX.out.csi)

            def haplotypecaller_vcf_tbi = haplotypecaller_vcf
                .join(haplotypecaller_indices, failOnDuplicate:true, failOnMismatch: true)

            def final_vcf = Channel.empty()

            //
            // MODULE: VariantFiltration from GATK4
            // Filter variant calls based on certain criteria
            //
            if (!params.skip_variantfiltration && !params.bam_csi_index ) {

                GATK4_VARIANTFILTRATION(
                    haplotypecaller_vcf_tbi,
                    fasta,
                    fasta_fai,
                    dict
                )

                def filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
                final_vcf    = filtered_vcf
                ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)
            } else {
                final_vcf = haplotypecaller_vcf
            }

            //
            // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
            //
            if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff') || params.annotate_tools.contains('vep'))) {

                def vep_fasta = [[id: 'null'], []]

                if (params.vep_include_fasta) {
                    vep_fasta = fasta
                }

                VCF_ANNOTATE_ALL(
                    final_vcf.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                    vep_fasta,
                    params.annotate_tools,
                    params.snpeff_genome ? "${params.snpeff_genome}.${params.snpeff_db}" : "${params.genome}.${params.snpeff_db}",
                    snpeff_cache,
                    vep_genome,
                    vep_species,
                    vep_cache_version,
                    vep_cache,
                    vep_extra_files,
                )

                // Gather QC ch_reports
                ch_reports  = ch_reports.mix(VCF_ANNOTATE_ALL.out.reports)
                ch_versions = ch_versions.mix(VCF_ANNOTATE_ALL.out.versions)
            }

        } else {

            //
            // MODULE: CombineGVCFS from GATK4
            // Merge multiple GVCF files into one GVCF
            //
            GATK4_COMBINEGVCFS(
                haplotypecaller_out,
                fasta.map { _meta, fasta_ -> fasta_ },
                fasta_fai.map { _meta, fai -> fai },
                dict.map { _meta, dict_ -> dict_ }
            )
            def haplotypecaller_gvcf = GATK4_COMBINEGVCFS.out.combined_gvcf
            ch_versions  = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions)

            //
            // MODULE: Index the VCF using TABIX
            //
            TABIXGVCF(haplotypecaller_gvcf)

            ch_versions  = ch_versions.mix(TABIXGVCF.out.versions)

        }
    }

    //
    // Collate and save software ch_versions
    //
    def collated_versions = softwareVersionsToYAML(ch_versions).collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    //
    def val_multiqc_report = Channel.empty()

    if (!params.skip_multiqc){
        def multiqc_files =  Channel.empty()

        def multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        def multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        def multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)   : Channel.empty()
        def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        def multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        def methods_description                = Channel.value(methodsDescriptionText(multiqc_custom_methods_description))

        multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(collated_versions)
        multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC (
            multiqc_files.collect(),
            multiqc_config.toList(),
            multiqc_custom_config.toList(),
            multiqc_logo.toList(),
            [],
            []
        )
        val_multiqc_report = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    multiqc_report = val_multiqc_report     // channel: /path/to/multiqc_report.html
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
