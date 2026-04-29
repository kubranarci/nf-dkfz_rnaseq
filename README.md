# nf/dkfz_rnaseq

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf/dkfz_rnaseq)
[![GitHub Actions CI Status](https://github.com/nf/dkfz_rnaseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf/dkfz_rnaseq/actions/workflows/nf-test.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

**nf/dkfz_rnaseq** is a bioinformatics pipeline that processes bulk RNA sequencing data (both single-end and paired-end) to perform alignment, quality control, gene expression quantification, and fusion transcript detection. It ingests FASTQ files and standard genomic references (FASTA, GTF) to dynamically map reads using STAR, quantify transcript abundances using featureCounts and Kallisto, and detect structural fusions using Arriba. The pipeline accurately routes data through a highly parallelized Nextflow architecture to output sorted alignments, extensive QC metrics, raw count matrices, rescaled TPMs, and annotated fusion calls.

### Pipeline Steps

1. **Index Generation**: Dynamic building of indices for mapping (`STAR`, `Kallisto`)
2. **Alignment**: 2-pass read mapping and chimeric alignment generation (`STAR`)
3. **Post-Alignment Processing**: Duplicate marking, sorting, indexing, and flagstat calculation (`Sambamba` and `Samtools`)
4. **Quality Control**: Comprehensive metric calculation (`RNA-SeQC`, `QualiMap2`)
5. **Gene/Exon Counting**: Read counting for standard gene models and exonic parts (`featureCounts`)
6. **Transcript Quantification**: Pseudoalignment, abundance quantification, and TPM rescaling (`Kallisto`)
7. **Fusion Detection**: Fusion transcript discovery and structural visualization (`Arriba`)
8. **QC Aggregation**: Consolidated JSON report generation from Sambamba and RNA-SeQC metrics (`Custom Perl`)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that strictly adheres to the pipeline's JSON schema. 

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,status,lane,replica,library,platform_unit
patient_01,patient_01_tumor_R1.fastq.gz,patient_01_tumor_R2.fastq.gz,tumor,L001,1,libA,FC123_L1
patient_01,patient_01_ctrl_R1.fastq.gz,patient_01_ctrl_R2.fastq.gz,control,L001,1,libB,FC123_L1
patient_02,patient_02_tumor_R1.fastq.gz,,tumor,L002,1,libC,FC124_L2
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired-end). For single-end data (as shown in the third row), leave the fastq_2 column empty. The pipeline uses the sample, status, and replica columns to automatically group multiple lanes belonging to the same read group before alignment.

Now, you can run the pipeline using the following minimal command:

```bash
nextflow run nf/dkfz_rnaseq \
   -profile <docker/singularity/conda> \
   --input samplesheet.csv \
   --fasta <path_to_reference.fa> \
   --gtf <path_to_annotation.gtf> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits


nf/dkfz_rnaseq was ported to Nextflow by kubranarci.

This pipeline is a modern Nextflow adaptation of the original HIPO2 RNAseq workflow developed at the German Cancer Research Center (DKFZ). We thank the original authors for their extensive assistance and foundational architecture:

- Naveed Ishaque

- Michael Heinold

- Jeongbin Park

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf/dkfz_rnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
