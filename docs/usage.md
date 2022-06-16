# nf-core/rnavar: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnavar/usage](https://nf-co.re/rnavar/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## RUNNING THE PIPELINE

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/rnavar --input samplesheet.csv --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/rnavar
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnavar releases page](https://github.com/nf-core/rnavar/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,unstranded
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,unstranded
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,unstranded
```

#### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,reverse
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,,reverse
```

| Column         | Description                                                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandedness` | Sample strand-specificity. Must be one of `unstranded`, `forward` or `reverse`.                                                                                                        |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## PIPELINE PARAMETERS AND DESCRIPTION

## Reference genome files

The minimum reference genome requirements are a FASTA and GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally.

> **NB:** Compressed reference files are also supported by the pipeline i.e. standard files with the `.gz` extension and indices folders with the `tar.gz` extension.

The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--star_index '/path/to/STAR/index/'`) or via a custom config file.

> **NB:** If you are supplying a pre-built genome index file via `--star_index`, please ensure that the index has been generated with the latest STAR version i.e. v2.7.9a or above. In case if the pipeline found an incompatible index, it will generate a new one using the reference genome which will consume time and memory unnecessarily.

- If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.
- If `--gff` is provided as input then this will be converted to a GTF file, or the latter will be used if both are provided.
- The `--exon_bed` parameter file is expected to be exon coordinates with at least three columns i.e., <chrom> <exon_position_start> <exon_position_end> in the file. The <exon_postion_start> should be 0-based. If this parameter is not provided, the exon coordinates are extracted from the GTF file and generates a bed file by the process `GTF2BED`.
- If `--star_index` is not provided then it will be generated from the reference genome FASTA file using `STAR --runmode genomeGenerate` command.

> **NB:** In case if you are providing a GTF and/or a BED file, please ensure that the chromosomes and contigs in the files are also present in the genome FASTA (and in the .dict) file. Otherwise `GATK BedToIntervalList` module is likely to fail if the chromosomes/contigs do not match with the reference genome data.

### Recommendation when using very large genomes

When the pipeline is used on very large genomes having chromosome size greater than 512Mb (e.g. Chromosome 1 of Monodelphis domestica has a size of 748055161bp), please make sure that `--bam_csi_index` parameter is provided in order to use coordinate sorted index (CSI) instead of standard binary alignment index (BAI).

> **NB:** When `--bam_csi_index` is used, variant filtration step will be disabled as `GATK VariantFiltration` does not currently support CSI index for the input VCF. It may be incorporated in the future when newer GATK versions support CSI for VCF inputs.

## Alignment options

The pipeline uses [STAR](https://github.com/alexdobin/STAR) to map the raw FastQ reads to the reference genome. STAR is fast but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

By default, STAR runs in `2-pass` mode. For the most sensitive novel junction discovery, it is recommend running STAR in the 2-pass
mode. It does not increase the number of detected novel junctions, but allows to detect more splices reads mapping to novel junctions. The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the first pass, and use them as ”annotated” junctions for the 2nd pass mapping. You can turn off this feature by setting `--star_twopass false` in command line.

Read length is an important parameter therefore it has to be used carefully. The default is set to 150, but it has to be changed according to the input reads. For example, if the input read length is 2x151bp, then you use `--read_length 151`. The `--read_length` parameter is used while generating an index as well as in the alignment process. In both processes, the pipeline use (read_length - 1) to the STAR parameter `--sjdbOverhang` as recommended in STAR documentation.

> **NB:** Read length `--read_length` is an important parameter, therefore it has to set according to the input read length. If you are supplying a pre-built genome index, please make sure that you have used the same (read_length -1) during the genomeGenerate step.

STAR alignment generates a coordinated-sorted BAM file as output. The coordinate-sorting process can be very memory intensive when the input data is deep sequenced or the genome has many highly expressed loci. When the pipeline runs on memory constrained environment, sorting step may fail due to low memory. In such cases you may adjust the limit parameters such as `--star_limitBAMsortRAM`, `--star_outBAMsortingBinsN` and `--star_limitOutSJcollapsed` to increase the sorting memory and genomic bins. Refer the parameter documentation for the default values and adjust as appropriate based on your memory availability.

## Preprocessing options

Marking duplicate reads is performed using `GATK4 MarkDuplicates` tool. The tool does not remove duplicate reads by default, however you can set `--remove_duplicates true` to remove them.

GATK best practices has been followed in this pipeline for RNA analysis, hence it uses GATK modules such as `SplitNCigarReads`, `BaseRecalibrator`, `ApplyBQSR`. The `BaseRecalibrator` process requires known variants sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation.You can supply the VCF and index files using parameters such as `--dbsnp`, `--dbsnp_tbi`, `--known_indels`, `--known_indels_tbi`.

> **NB:** Base recalibration can be turned off using `--skip_baserecalibration true` option. This is useful when you are analyzing data from non-model organisms where there is no known variant datasets exist.

`GATK SplitNCigarReads` is very time consuming step, therefore we made an attempt to break the GTF file into multiple chunks (scatters) using `GATK IntervalListTools` to run the process independently on each chunk in a parallel way to speed up the analysis. The default number of splits is set to 25, that means the GTF file is split into 25 smaller files and run `GATK SplitNCigarReads` on each of them in parallel. You can modify the number of splits using parameter `--gatk_interval_scatter_count`.

## Variant calling and filtering

`GATK HaplotypeCaller` is used for variant calling with default minimum phred-scaled confidence threshold as 20. This value can be changed using paramerter `--gatk_hc_call_conf`.

The pipeline runs a hard-filtering step on the variants by default. It does not filter out any variants, rather it flags i.e. PASS or other flags such as FS, QD, SnpCluster, etc. in FILTER column of the VCF. The following are the default filter criteria, however it can be changed using the respective parameters.

- `--gatk_vf_cluster_size` is set to 3. It is the number of SNPs which make up a cluster.
- `--gatk_vf_window_size` is set to 35. The window size (in bases) in which to evaluate clustered SNPs.
- `--gatk_vf_fs_filter` is set to 30.0. Filter based on FisherStrand > 30.0. It is the Phred-scaled probability that there is strand bias at the site.
- `--gatk_vf_qd_filter` is set to 2.0 meaning filter variants if Quality By Depth filter is < 2.0.

Variant filtering is an optional step. You can skip it using `--skip_variantfiltration` parameter.

## Variant annotation

The annotation of variants is performed using snpEff and VEP. The parameter to use is `--annotate_tools snpeff` or `--annotate_tools vep`. You can even run both snpEff and VEP using `--annotate_tools merge`, in this case the output VCF file will have both snpEff and VEP annotations combined.

You can skip the variant annotation step using `--skip_variantannotation` parameter or without passing `--annotate_tools` options.

### Annotation cache

Both `snpEff` and `VEP` enable usage of cache.
If cache is available on the machine where `rnavar` is run, it is possible to run annotation using cache.
You need to specify the cache directory using `--snpeff_cache` and `--vep_cache` in the command lines or within configuration files.
The cache will only be used when `--annotation_cache` and cache directories are specified (either in command lines or in a configuration file).

Example:

```bash
nextflow run nf-core/rnavar --input samplesheet.csv --genome GRCh38 -profile docker --annoate_tools snpEff --snpeff_cache </path/to/snpEff/cache> --annotation_cache
nextflow run nf-core/rnavar --input samplesheet.csv --genome GRCh38 -profile docker --annotate_tools VEP --vep_cache </path/to/VEP/cache> --annotation_cache
```

### Download annotation cache

A `Nextflow` helper script [link](https://raw.githubusercontent.com/nf-core/sarek/master/download_cache.nf) has been designed to help downloading `snpEff` and `VEP` caches.
Such files are meant to be shared between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --snpeff_cache </path/to/snpEff/cache> --snpeff_db <snpEff DB version> --genome <GENOME>
nextflow run download_cache.nf --vep_cache </path/to/VEP/cache> --species <species> --vep_cache_version <VEP cache version> --genome <GENOME>
```

### Using VEP CADD plugin

To enable the use of the `VEP` `CADD` plugin:

- Download the `CADD` files
- Specify them (either on the command line, like in the example or in a configuration file)
- use the `--cadd_cache` flag

Example:

```bash
nextflow run nf-core/rnavar --input samplesheet.csv --genome GRCh38 -profile docker --annotate_tools VEP VEP --cadd_cache \
    --cadd_indels </path/to/CADD/cache/InDels.tsv.gz> \
    --cadd_indels_tbi </path/to/CADD/cache/InDels.tsv.gz.tbi> \
    --cadd_wg_snvs </path/to/CADD/cache/whole_genome_SNVs.tsv.gz> \
    --cadd_wg_snvs_tbi </path/to/CADD/cache/whole_genome_SNVs.tsv.gz.tbi>
```

### Downloading CADD files

An helper script has been designed to help downloading `CADD` files.
Such files are meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --cadd_cache </path/to/CADD/cache> --cadd_version <CADD version> --genome <GENOME>
```

## GENERAL NEXTFLOW ARGUMENTS

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## CUSTOM CONFIGURATIONS

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnavar/blob/master/conf/base.config#L17) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnavar pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNAVAR:RNAVAR:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNAVAR:RNAVAR:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNAVAR:RNAVAR:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process.
The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnavar Github repo](https://github.com/nf-core/rnavar/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnavar/blob/master/modules/nf-core/modules/star/align/main.nf#L3).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnavar/blob/master/conf/base.config#L37-L41) which in this case is defined as `72GB`.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least `72GB` of memory, in this case increased to `100GB`.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNAVAR:RNAVAR:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNAVAR:RNAVAR:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
