# RNAseq-vax
The aim of this pipeline is to processed bam files from NGI-RNAseq and output variants and allele specific expression

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* add read groups ( via picard)
* splitCigarReads (via GATK)
* call variants ( via GATK HaplotypeCaller)
* Variant filtering (via GATK)
* annotate variants for further downstream analysis (via VEP)
* select variants (via GATK)
* allele specific expression (via GATK)

## output

Directories output
* haplotypeCaller (raw variants, from haplotypecaller, bziped and indexed)
* VariantFiltration (hard filters applied according to https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq )
* BiallelecVCF (only bialleleic SNVs, for input into allele specific expression)
* VEP( annotated vcf files)
* AlleleSpecificExpression ( allele specific expression for each sample)
