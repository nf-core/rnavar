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
    ch_input
    ch_dbsnp
    ch_dbsnp_tbi
    ch_dict
    ch_exon_bed
    ch_fasta
    ch_fasta_fai
    ch_gtf
    ch_known_indels
    ch_known_indels_tbi
    ch_star_index
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

    // To gather used softwares versions for MultiQC
    def ch_versions = Channel.empty()

    // Parse the input data
    ch_parsed_input = ch_input.groupTuple().map{ samplesheet -> checkSamplesAfterGrouping(samplesheet) }
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
        ch_parsed_input.cram,
        ch_parsed_input.bam,
    )
    ch_versions = ch_versions.mix(PREPARE_ALIGNMENT.out.versions)

    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ(ch_parsed_input.multiple)

    def ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_parsed_input.single)

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // MODULE: Generate QC summary using FastQC
    FASTQC(ch_cat_fastq)
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{ _meta, logs -> logs })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
    //

    GATK4_BEDTOINTERVALLIST(ch_exon_bed, ch_dict)
    def ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    //
    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    //
    def ch_interval_list_split = Channel.empty()
    if (!params.skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(ch_interval_list)
        ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ _meta, bed -> [bed] }.flatten()
    }
    else {
        ch_interval_list_split = ch_interval_list.map { _meta, bed -> bed }
    }

    //
    // SUBWORKFLOW: Perform read alignment using STAR aligner
    //

    if (params.aligner == 'star') {
        FASTQ_ALIGN_STAR(ch_cat_fastq,
            ch_star_index,
            ch_gtf,
            params.star_ignore_sjdbgtf,
            seq_platform,
            seq_center,
            ch_fasta,
            [[:],[]]) //ch_transcripts_fasta)

        def ch_genome_bam    = FASTQ_ALIGN_STAR.out.bam

        // Gather QC reports
        ch_reports           = ch_reports.mix(FASTQ_ALIGN_STAR.out.log_out)
        ch_reports           = ch_reports.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_versions          = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Mark duplicates with GATK4
        //
        BAM_MARKDUPLICATES_PICARD(ch_genome_bam,
            ch_fasta,
            ch_fasta_fai)

        def ch_genome_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(BAM_MARKDUPLICATES_PICARD.out.bai, remainder: true)
            .join(BAM_MARKDUPLICATES_PICARD.out.csi, remainder: true) // TODO fix this bottleneck
            .map{meta, bam, bai, csi ->
                if (bai) [meta, bam, bai]
                else [meta, bam, csi]
            }
            .mix(PREPARE_ALIGNMENT.out.bam)

        //Gather QC reports
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]))
        ch_reports                = ch_reports.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]))
        ch_versions               = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        //
        // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string(e.g. spanning splicing events in RNAseq data).
        //

        SPLITNCIGAR(ch_genome_bam_bai,
            ch_fasta,
            ch_fasta_fai.map { _meta, fai -> fai },
            ch_dict,
            ch_interval_list_split
        )

        def ch_splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
        ch_versions                 = ch_versions.mix(SPLITNCIGAR.out.versions)

        //
        // MODULE: BaseRecalibrator from GATK4
        // Generates a recalibration table based on various co-variates
        //
        def ch_bam_variant_calling = Channel.empty()

        if (!params.skip_baserecalibration) {
            // known_sites is made by grouping both the dbsnp and the known indels ressources
            // they can either or both be optional
            def ch_known_sites     = ch_dbsnp.concat(ch_known_indels).collect()
            def ch_known_sites_tbi = ch_dbsnp_tbi.concat(ch_known_indels_tbi).collect()

            def ch_interval_list_recalib = ch_interval_list.map{ _meta, bed -> [bed] }.flatten()
            def ch_splitncigar_bam_bai_interval = ch_splitncigar_bam_bai.combine(ch_interval_list_recalib)
                .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval] }

            GATK4_BASERECALIBRATOR(
                ch_splitncigar_bam_bai_interval,
                ch_fasta.map{ _meta, fasta -> [fasta] },
                ch_fasta_fai.map { _meta, fai -> fai },
                ch_dict.map{ _meta, dict -> [dict] },
                ch_known_sites,
                ch_known_sites_tbi
            )
            def ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table

            // Gather QC reports
            ch_reports  = ch_reports.mix(ch_bqsr_table.map{ _meta, table -> table})
            ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

            def ch_bam_applybqsr       = ch_splitncigar_bam_bai.join(ch_bqsr_table)

            def ch_interval_list_applybqsr = ch_interval_list.map{ _meta, bed -> [bed] }.flatten()
            def ch_applybqsr_bam_bai_interval = ch_bam_applybqsr.combine(ch_interval_list_applybqsr)
                .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]}

            //
            // MODULE: ApplyBaseRecalibrator from GATK4
            // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
            //
            RECALIBRATE(
                params.skip_multiqc,
                ch_applybqsr_bam_bai_interval,
                ch_dict.map{ _meta, dict -> [dict] },
                ch_fasta_fai.map { _meta, fai -> fai },
                ch_fasta.map{ _meta, fasta -> [fasta] }
            )

            ch_bam_variant_calling = RECALIBRATE.out.bam

            // Gather QC reports
            ch_reports  = ch_reports.mix(RECALIBRATE.out.qc.collect{it[1]}.ifEmpty([]))
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions)
        } else {
            ch_bam_variant_calling = ch_splitncigar_bam_bai
        }

        // Run haplotyper even in the absence of dbSNP files
        def ch_dbsnp_for_haplotypecaller = Channel.empty()
        def ch_dbsnp_for_haplotypecaller_tbi = Channel.empty()
        if (!params.dbsnp){
            ch_dbsnp_for_haplotypecaller = [[id:'null'], []]
            ch_dbsnp_for_haplotypecaller_tbi = [[id:'null'], []]
        } else {
            ch_dbsnp_for_haplotypecaller     = ch_dbsnp.map{ vcf -> [[id:'dbsnp'], vcf] }
            ch_dbsnp_for_haplotypecaller_tbi = ch_dbsnp_tbi.map{ tbi -> [[id:'dbsnp'], tbi] }
        }

        ch_haplotypecaller_interval_bam = ch_bam_variant_calling.combine(ch_interval_list_split)
            .map{ meta, bam, bai, interval_list ->
                [ meta + [ id:meta.id + "_" + interval_list.baseName, sample:meta.id, variantcaller:'haplotypecaller' ], bam, bai, interval_list, [] ]
            }

        //
        // MODULE: HaplotypeCaller from GATK4
        // Calls germline SNPs and indels via local re-assembly of haplotypes.
        //

        GATK4_HAPLOTYPECALLER(
            ch_haplotypecaller_interval_bam,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_dbsnp_for_haplotypecaller,
            ch_dbsnp_for_haplotypecaller_tbi
        )

        def ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf.map{ meta, vcf -> [ meta + [id:meta.sample] - meta.subMap('sample'), vcf ] }.groupTuple() // TODO fix this bottleneck

        ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

        def ch_haplotypecaller_vcf = Channel.empty()
        if (!params.generate_gvcf){
            //
            // MODULE: MergeVCFS from GATK4
            // Merge multiple VCF files into one VCF
            //
            GATK4_MERGEVCFS(
                ch_haplotypecaller_raw,
                ch_dict
            )
            ch_haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
            ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

            //
            // MODULE: Index the VCF using TABIX
            //
            TABIX(
                ch_haplotypecaller_vcf
            )

            def ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_vcf
                .join(TABIX.out.tbi, by: [0], remainder: true)
                .join(TABIX.out.csi, by: [0], remainder: true) // TODO fix this bottleneck
                .map{meta, vcf, tbi, csi ->
                    if (tbi) [meta, vcf, tbi]
                    else [meta, vcf, csi]
                }

            ch_versions      = ch_versions.mix(TABIX.out.versions)
            def ch_final_vcf = Channel.empty()

            //
            // MODULE: VariantFiltration from GATK4
            // Filter variant calls based on certain criteria
            //
            if (!params.skip_variantfiltration && !params.bam_csi_index ) {

                GATK4_VARIANTFILTRATION(
                    ch_haplotypecaller_vcf_tbi,
                    ch_fasta,
                    ch_fasta_fai,
                    ch_dict
                )

                def ch_filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
                ch_final_vcf    = ch_filtered_vcf
                ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)
            } else {
                ch_final_vcf = ch_haplotypecaller_vcf
            }

            //
            // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
            //
            if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff') || params.annotate_tools.contains('vep'))) {

                def ch_vep_fasta = [[id: 'null'], []]

                if (params.vep_include_fasta) {
                    ch_vep_fasta = ch_fasta
                }

                VCF_ANNOTATE_ALL(
                    ch_final_vcf.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                    ch_vep_fasta,
                    params.annotate_tools,
                    params.snpeff_genome ? "${params.snpeff_genome}.${params.snpeff_db}" : "${params.genome}.${params.snpeff_db}",
                    snpeff_cache,
                    vep_genome,
                    vep_species,
                    vep_cache_version,
                    vep_cache,
                    vep_extra_files,
                )

                // Gather QC reports
                ch_reports  = ch_reports.mix(VCF_ANNOTATE_ALL.out.reports)
                ch_versions = ch_versions.mix(VCF_ANNOTATE_ALL.out.versions)
            }

        } else {
            def ch_combinegvcfs_input = GATK4_HAPLOTYPECALLER.out.vcf
                .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch:true, failOnDuplicate:true)
                .map{ meta, vcf, tbi ->
                    def new_meta = meta + [id:meta.sample]
                    [new_meta, vcf, tbi]
                }
                .groupTuple() // TODO fix this bottleneck

            //
            // MODULE: CombineGVCFS from GATK4
            // Merge multiple GVCF files into one GVCF
            //
            GATK4_COMBINEGVCFS(
                ch_combinegvcfs_input,
                ch_fasta.map { _meta, fasta -> fasta },
                ch_fasta_fai.map { _meta, fai -> fai },
                ch_dict.map { _meta, dict -> dict }
            )
            def ch_haplotypecaller_gvcf = GATK4_COMBINEGVCFS.out.combined_gvcf
            ch_versions  = ch_versions.mix(GATK4_COMBINEGVCFS.out.versions)

            //
            // MODULE: Index the VCF using TABIX
            //
            TABIXGVCF(ch_haplotypecaller_gvcf)

            ch_versions  = ch_versions.mix(TABIXGVCF.out.versions)

        }
    }

    //
    // Collate and save software versions
    //
    def ch_collated_versions = softwareVersionsToYAML(ch_versions).collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    //
    def val_multiqc_report = Channel.empty()

    if (!params.skip_multiqc){
        def ch_multiqc_files =  Channel.empty()

        def ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        def ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        def ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)   : Channel.empty()
        def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
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
