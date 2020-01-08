# RNAseq-vax
The aim of this pipeline is to processed bam files from nf-core/RNAseq and output variants and allele specific expression

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
Following steps are performed:

* add read groups ( via picard)
* splitCigarReads (via GATK)
* call variants ( via GATK HaplotypeCaller)
* Variant filtering (via GATK)
* select variants (via GATK)
* allele specific expression (via GATK)


## Output folders and files

Directories output
* HaplotypeCaller (raw variants, from haplotypecaller, bziped and indexed)
* VariantFiltration (hard filters applied according to [GATK Best Practices workflow for SNP and indel calling on RNAseq data](https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq)
* BiallelecVCF (only bialleleic SNVs, for input into allele specific expression)

*Please note that these files are not filtered and that needs to be done as well*