# RNA-seq_AVC
this pipeline takes already aligned data from NGI-RNAseq pipeline and calls variants with
haplotypecaller as wella as allele specific variants

[![Build Status](https://travis-ci.org/rna-seq_AVC.svg?branch=master)](https://travis-ci.org/rna-seq_AVC)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)


### Introduction
RNA-seq_AVC: this pipeline takes already aligned data from NGI-RNAseq pipeline and calls variants and allele specifice variant expressions

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The RNA-seq_AVC pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits
This pipeline was written by Aron T. Skaftason ([arontommi](https://github.com/arontommi)) but mostly cannibalised from [NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq/)

