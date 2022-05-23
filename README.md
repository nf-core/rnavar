# ![nf-core/rnavar](docs/images/nf-core-rnavar_logo_light.png#gh-light-mode-only) ![nf-core-rnavar](docs/images/nf-core/rnavar_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/rnavar/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rnavar/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnavar/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rnavar/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/rnavar/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/rnavar)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnavar-4A154B?logo=slack)](https://nfcore.slack.com/channels/rnavar)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnavar** is a bioinformatics best-practice analysis pipeline for GATK4 RNA variant calling.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/rnavar/results).

## Pipeline summary

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align reads to reference genome ([`STAR`](https://github.com/alexdobin/STAR))
4. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Duplicate read marking ([`GATK4 MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard))
6. Splits reads that contain Ns in their cigar string ([`GATK4 SplitNCigarReads`](https://gatk.broadinstitute.org/hc/en-us/articles/4409917482651-SplitNCigarReads))
7. Estimate and correct systematic bias using base quality score recalibration ([`GATK4 BaseRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/4409897206043-BaseRecalibrator), [`GATK4 ApplyBQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/4409897168667-ApplyBQSR))
8. Convert a BED file to a Picard Interval List ([`GATK4 BedToIntervalList`](https://gatk.broadinstitute.org/hc/en-us/articles/4409924780827-BedToIntervalList-Picard-))
9. Scatter one interval-list into many interval-files ([`GATK4 IntervalListTools`](https://gatk.broadinstitute.org/hc/en-us/articles/4409917392155-IntervalListTools-Picard-))
10. Call SNPs and indels ([`GATK4 HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/4409897180827-HaplotypeCaller))
11. Merge multiple VCF files into one VCF ([`GATK4 MergeVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/4409924817691-MergeVcfs-Picard-))
12. Index the VCF ([`Tabix`](http://www.htslib.org/doc/tabix.html))
13. Filter variant calls based on certain criteria ([`GATK4 VariantFiltration`](https://gatk.broadinstitute.org/hc/en-us/articles/4409897204763-VariantFiltration))
14. Annotate variants ([`snpEff`](https://pcingola.github.io/SnpEff/se_introduction/), [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html))
15. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

### Summary of tools and version used in the pipeline:

| Tool        | Version |
| ----------- | ------- |
| FastQC      | 0.11.9  |
| STAR        | 2.7.9a  |
| Samtools    | 1.15.1  |
| GATK        | 4.2.6.1 |
| Tabix       | 1.11    |
| SnpEff      | 5.0     |
| Ensembl VEP | 104.3   |
| MultiQC     | 1.12    |

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/rnavar -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run nf-core/rnavar -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --genome GRCh38
   ```

## Documentation

The nf-core/rnavar pipeline comes with documentation about the pipeline [usage](https://nf-co.re/rnavar/usage), [parameters](https://nf-co.re/rnavar/parameters) and [output](https://nf-co.re/rnavar/output).

## Credits

These scripts were originally written in Nextflow DSL2 for use at the [Barntumörbanken, Karolinska Institutet](https://ki.se/forskning/barntumorbanken), by Praveen Raj ([@praveenraj2018](https://github.com/praveenraj2018)) and Maxime U. Garcia ([@maxulysse](https://github.com/maxulysse)).

The pipeline is primarily maintained by Praveen Raj ([@praveenraj2018](https://github.com/praveenraj2018)) and Maxime U. Garcia ([@maxulysse](https://github.com/maxulysse)) from [Barntumörbanken, Karolinska Institutet](https://ki.se/forskning/barntumorbanken).

Many thanks to other who have helped out along the way too, including (but not limited to):
[@ewels](https://github.com/ewels),
[@drpatelh](https://github.com/drpatelh).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnavar` channel](https://nfcore.slack.com/channels/rnavar) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
