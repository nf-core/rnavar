/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core
include { CAT_FASTQ                 } from '../modules/nf-core/cat/fastq'
include { FASTQC                    } from '../modules/nf-core/fastqc'
include { GATK4_BASERECALIBRATOR    } from '../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_BEDTOINTERVALLIST   } from '../modules/nf-core/gatk4/bedtointervallist'
include { GATK4_COMBINEGVCFS        } from '../modules/nf-core/gatk4/combinegvcfs'
include { GATK4_HAPLOTYPECALLER     } from '../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_INTERVALLISTTOOLS   } from '../modules/nf-core/gatk4/intervallisttools'
include { GATK4_MERGEVCFS           } from '../modules/nf-core/gatk4/mergevcfs'
include { GATK4_VARIANTFILTRATION   } from '../modules/nf-core/gatk4/variantfiltration'
include { MOSDEPTH                  } from '../modules/nf-core/mosdepth'
include { SEQ2HLA                   } from '../modules/nf-core/seq2hla'
include { TABIX_TABIX as TABIX      } from '../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIXGVCF  } from '../modules/nf-core/tabix/tabix'
include { UMITOOLS_EXTRACT          } from '../modules/nf-core/umitools/extract'

// local
include { PREPARE_ALIGNMENT         } from '../subworkflows/local/prepare_alignment'
include { RECALIBRATE               } from '../subworkflows/local/recalibrate'
include { SPLITNCIGAR               } from '../subworkflows/local/splitncigar'
include { VCF_ANNOTATE_ALL          } from '../subworkflows/local/vcf_annotate_all'

// nf-core
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { FASTQ_ALIGN_STAR          } from '../subworkflows/nf-core/fastq_align_star'

// local
include { checkSamplesAfterGrouping } from '../subworkflows/local/utils_nfcore_rnavar_pipeline'

/*
========================================================================================
    RUN MAIN WORKFLOW RNAVAR
========================================================================================
*/

workflow RNAVAR {
    take:
    input
    bcftools_annotations
    bcftools_annotations_tbi
    bcftools_columns
    bcftools_header_lines
    dbsnp
    dbsnp_tbi
    dict
    exon_bed
    fasta
    fasta_fai
    gtf
    known_sites
    known_sites_tbi
    star_index
    snpeff_cache
    snpeff_db
    vep_genome
    vep_species
    vep_cache_version
    vep_include_fasta
    vep_cache
    vep_extra_files
    aligner
    star_ignore_sjdbgtf
    tools

    main:
    // Parse the input data
    parsed_input = input
        .groupTuple()
        .map { samplesheet -> checkSamplesAfterGrouping(samplesheet) }
        .branch { meta, fastqs, bam, bai, cram, crai, vcf, tbi ->
            single: fastqs.size() == 1
            return [meta, fastqs.flatten()]
            multiple: fastqs.size() > 1
            return [meta, fastqs.flatten()]
            bam: bam
            return [meta, bam, bai]
            cram: cram
            return [meta, cram, crai]
            vcf: vcf
            return [meta, vcf, tbi]
        }

    // MODULE: Prepare the alignment files (index BAM/CRAM files that are missing an index)
    PREPARE_ALIGNMENT(parsed_input.bam, parsed_input.cram)

    MOSDEPTH(parsed_input.cram.map { meta, cram, crai -> [meta, cram, crai, []] }, fasta, false)

    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ(parsed_input.multiple)

    def reads_input_all = CAT_FASTQ.out.reads.mix(parsed_input.single)

    // MODULE: Generate QC summary using FastQC
    FASTQC(reads_input_all)

    // MODULE: Extract UMIs from reads
    UMITOOLS_EXTRACT(reads_input_all.filter { 'umitools' in tools })

    def reads_preprocessed = 'umitools' in tools ? UMITOOLS_EXTRACT.out.reads : reads_input_all

    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
    GATK4_BEDTOINTERVALLIST(exon_bed, dict)

    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list.filter { 'intervallisttools' in tools })

    def interval_list_split = 'intervallisttools' in tools
        ? GATK4_INTERVALLISTTOOLS.out.interval_list.map { _meta, bed -> [bed] }.collect()
        : GATK4_BEDTOINTERVALLIST.out.interval_list.map { _meta, bed -> [[bed]] }.collect()

    // MODULE: HLATyping with Seq2HLA
    SEQ2HLA(reads_preprocessed.filter { 'seq2hla' in tools })

    // SUBWORKFLOW: Perform read alignment using STAR aligner

    if (aligner == 'star') {
        FASTQ_ALIGN_STAR(
            reads_preprocessed,
            star_index,
            gtf,
            star_ignore_sjdbgtf,
            fasta.join(fasta_fai).collect(),
            [[:], []],
        )

        // SUBWORKFLOW: Mark duplicates with GATK4
        BAM_MARKDUPLICATES_PICARD(FASTQ_ALIGN_STAR.out.bam, fasta.join(fasta_fai).collect())

        def genome_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(BAM_MARKDUPLICATES_PICARD.out.index, failOnDuplicate: true, failOnMismatch: true)
            .mix(PREPARE_ALIGNMENT.out.reads_index)

        // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
        // Splits reads that contain Ns in their cigar string(e.g. spanning splicing events in RNAseq data).

        SPLITNCIGAR(
            genome_bam_bai,
            fasta,
            fasta_fai,
            dict,
            interval_list_split,
        )

        // MODULE: BaseRecalibrator from GATK4
        // Generates a recalibration table based on various co-variates
        def bam_variant_calling = channel.empty()

        if ('baserecalibrator' in tools) {
            def splitncigar_bam_bai_interval = SPLITNCIGAR.out.bam_bai.combine(GATK4_BEDTOINTERVALLIST.out.interval_list.map { _meta, bed -> [bed] }.flatten())

            GATK4_BASERECALIBRATOR(
                splitncigar_bam_bai_interval,
                fasta,
                fasta_fai,
                dict,
                known_sites,
                known_sites_tbi,
            )

            def bam_applybqsr = SPLITNCIGAR.out.bam_bai.join(GATK4_BASERECALIBRATOR.out.table)

            def applybqsr_bam_bai_interval = bam_applybqsr
                .combine(GATK4_BEDTOINTERVALLIST.out.interval_list.map { _meta, bed -> [bed] }.flatten())
                .map { meta, bam, bai, table, interval -> [meta, bam, bai, table, interval] }

            // MODULE: ApplyBaseRecalibrator from GATK4
            // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
            RECALIBRATE(
                applybqsr_bam_bai_interval,
                dict.map { _meta, dict_ -> [dict_] },
                fasta_fai,
                fasta,
            )

            bam_variant_calling = RECALIBRATE.out.bam
        }
        else {
            bam_variant_calling = SPLITNCIGAR.out.bam_bai
        }

        def haplotypecaller_interval_bam = bam_variant_calling
            .combine(interval_list_split)
            .map { meta, bam, bai, interval_lists -> [meta + [interval_count: interval_lists.size()], bam, bai, interval_lists] }
            .transpose(by: 3)
            .map { meta, bam, bai, interval_list_ -> [meta + [id: meta.id + "_" + interval_list_.baseName, sample: meta.id, variantcaller: 'haplotypecaller'], bam, bai, interval_list_, []] }

        // MODULE: HaplotypeCaller from GATK4
        // Calls germline SNPs and indels via local re-assembly of haplotypes.

        GATK4_HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
        )

        def haplotypecaller_out = GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi -> [groupKey(meta + [id: meta.sample] - meta.subMap('sample', "interval_count"), meta.interval_count), vcf, tbi] }
            .groupTuple()

        if ('combinegvcfs' in tools) {
            // MODULE: CombineGVCFS from GATK4
            // Merge multiple GVCF files into one GVCF
            GATK4_COMBINEGVCFS(
                haplotypecaller_out,
                fasta.map { _meta, fasta_ -> fasta_ },
                fasta_fai.map { _meta, fai -> fai },
                dict.map { _meta, dict_ -> dict_ },
            )

            // MODULE: Index the VCF using TABIX
            TABIXGVCF(GATK4_COMBINEGVCFS.out.combined_gvcf.map { meta, vcf -> [meta, vcf, [], []] })
        }
        else {
            // MODULE: MergeVCFS from GATK4
            // Merge multiple VCF files into one VCF
            def haplotypecaller_raw = haplotypecaller_out.map { meta, vcfs, _tbis -> [meta, vcfs] }
            GATK4_MERGEVCFS(
                haplotypecaller_raw,
                dict,
            )

            // MODULE: Index the VCF using TABIX
            TABIX(GATK4_MERGEVCFS.out.vcf.map { meta, vcf -> [meta, vcf, [], []] })

            def haplotypecaller_vcf_tbi = GATK4_MERGEVCFS.out.vcf.join(TABIX.out.index, failOnDuplicate: true, failOnMismatch: true)

            def final_vcf = channel.empty()

            // MODULE: VariantFiltration from GATK4
            // Filter variant calls based on certain criteria
            if ('variantfiltration' in tools) {
                GATK4_VARIANTFILTRATION(
                    haplotypecaller_vcf_tbi,
                    fasta,
                    fasta_fai,
                    dict,
                    [[:], []],
                )

                def filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
                final_vcf = filtered_vcf
            }
            else {
                final_vcf = GATK4_MERGEVCFS.out.vcf
            }

            // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
            if ('bcfann' in tools || 'merge' in tools || 'snpeff' in tools || 'vep' in tools) {

                final_vcf = final_vcf.mix(parsed_input.vcf.map { meta, vcf, _tbi -> [meta, vcf] })

                VCF_ANNOTATE_ALL(
                    final_vcf.map { meta, vcf -> [meta + [file_name: vcf.baseName], vcf] },
                    fasta.map { meta, fasta_ -> [meta, vep_include_fasta ? fasta_ : []] },
                    tools,
                    snpeff_db,
                    snpeff_cache,
                    vep_genome,
                    vep_species,
                    vep_cache_version,
                    vep_cache,
                    vep_extra_files,
                    bcftools_annotations,
                    bcftools_annotations_tbi,
                    bcftools_columns,
                    bcftools_header_lines,
                )
            }
        }
    }
}
