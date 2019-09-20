
please note : This pipeline is still in alpha stage. please use at your own risk


# RNAseq-VAX - RNAseq Variante Allele eXpression 

### Introduction

The aim of this pipeline ist to take of where [nf-core/rnaseq](https://github.com/nf-core/rnaseq) leaves and add additional informations on top. The pipeline follows [GATK guideline for variant calling in RNAseq](https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail). It then uses [ASEReadCounter](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) to retrive allele specific expression.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The RNAseq_VAX pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits
This pipeline was written by Aron T. Skaftason ([arontommi](https://github.com/arontommi)) but mostly cannibalised from [nf-core/rnaseq](https://github.com/nf-core/rnaseq)