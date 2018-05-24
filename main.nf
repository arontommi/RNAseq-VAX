#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         RNA-seq_AVC
========================================================================================
 RNA-seq_AVC Analysis Pipeline. Started 2018-02-06.
 #### Homepage / Documentation
 rna-seq_AVC
 #### Authors
 Aron T. Skaftason arontommi <aron.skaftason@ki.se> - https://github.com/arontommi>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     RNA-seq_AVC v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run RNAseq_AVC  -with-singularity rna-avc.img --project 'your_uppmax_project' --genome /sw/data/uppnex/reference/Homo_sapiens/hg19/program_files/GATK/concat.fasta
    
    Mandatory arguments:
      --genome                      Genome used to align 
      --project                     your Uppmax project
      -with-singularity             Singularity container


    """.stripIndent()
}



params.deduped_bam = 'T'
params.outdir = './results'
params.bamfolder = './results/markDuplicates/'
params.genome = false
params.project = false

params.rglb = '1'
params.rgpl = 'illumina'
params.rgpu = 'unit1'

wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")


if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if (params.bamfolder ) {
    Channel
        .fromPath(params.bamfolder+'*.bam')
        .set{bam_md}
}
else if ( !params.deduped_bam ){
    exit 1, "No Bam specified! do some aligning first!"
}

if ( params.genome ){
    if ( params.genome.endsWith('.fa')) {
        genomefasta = file(params.genome)
        genomefai = file(params.genome + '.fai')
        genomedict = file(params.genome - '.fa'+'.dict')
    }
    else if ( params.genome.endsWith('.fasta')) {
        genomefasta = file(params.genome)
        genomefai = file(params.genome + '.fai')
        genomedict = file(params.genome - '.fasta'+'.dict')
    }
}
/*
 * Readgroups added 
 */
process addReadGroups{ 
    tag "$bam_md.baseName"

    input:
    file bam_md
    val rglb from params.rglb
    val rgpl from params.rgpl
    val rgpu from params.rgpu

    output:
    set val("$name"), file("${name}.RG.bam"), file("${name}.RG.bam.bai") into rg_data

    script:
    if( task.memory == null ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    name = "${bam_md.baseName}"
    """
    java -Xmx${avail_mem}g -jar \$PICARD_HOME/picard.jar AddOrReplaceReadGroups \\
        I= $bam_md \\
        O= ${name}.RG.bam \\
        RGLB=$rglb \\
        RGPL=$rgpl \\
        RGPU=$rgpu \\
        RGSM=${bam_md.baseName}

    samtools index ${bam_md.baseName}.RG.bam
    """
}
/*
 * STEP 7 SplitNCigarReads
 */
process splitNCigarReads {
    tag "$rg_bam"

    input:
    set val(name), file(rg_bam), file(rg_bam_bai) from rg_data
    file genomefasta
    file genomefai
    file genomedict

    output:
    set val("$name"), file("${name}_split.bam"), file("${name}_split.bam.bai") into sc_data

    script:

    """
    java -jar \$GATK_HOME/gatk-package-4.0.1.2-local.jar SplitNCigarReads \\
    -R $genomefasta \\
    -I $rg_bam \\
    -O ${name}_split.bam 

    samtools index ${name}_split.bam
    """
}

/*
 * STEP 7 Haplotypecaller
 */
process haplotypeCaller {
    tag "$splitNCigar_bam.baseName"
    publishDir "${params.outdir}/haplotypeCaller", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".bam") || filename.endsWith(".bai")) null
                    else filename
                    }

    input:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai) from sc_data
    file genomefasta
    file genomefai
    file genomedict

    output:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file("${name}.vcf.gz"), file("${name}.vcf.gz.tbi") into ht_data

    script:

    """
    java -jar \$GATK_HOME/gatk-package-4.0.1.2-local.jar HaplotypeCaller \\
    -R $genomefasta \\
    -I $splitNCigar_bam \\
    --dont-use-soft-clipped-bases \\
    --standard-min-confidence-threshold-for-calling 20.0 \\
    -O ${name}.vcf 
    bgzip -c ${name}.vcf > ${name}.vcf.gz
    tabix -p vcf ${name}.vcf.gz
    """
}

/*
 * STEP 8 varfiltering
 */
process varfiltering {
    tag "$vcf.baseName"
    publishDir "${params.outdir}/VariantFiltration", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".bam") || filename.endsWith(".bai")) null
                    else filename
                    }


    input:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file(vcf), file(vcf_tbi) from ht_data
    file genomefasta
    file genomefai
    file genomedict


    output:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file("${name}.sorted.vcf.gz"), file("${name}.sorted.vcf.gz.tbi") into filtered_data
    
    script:

    """
    java -jar \$GATK_HOME/gatk-package-4.0.1.2-local.jar VariantFiltration \\
    -R $genomefasta \\
    -V $vcf \\
    -window 35 \\
    -cluster 3 \\
    -filter-name FS \\
    -filter "FS > 30.0" \\
    -filter-name QD \\
    -filter "QD < 2.0"  \\
    -O ${name}.sorted.vcf
    bgzip -c ${name}.sorted.vcf > ${name}.sorted.vcf.gz
    tabix -p vcf ${name}.sorted.vcf.gz
    """
}



 process selectvariants {
    tag "$filtered_vcf.baseName"
    publishDir "${params.outdir}/BiallelecVCF", mode: 'copy',
        saveAs: {filename ->
                    if (filename.endsWith(".bam") || filename.endsWith(".bai")) null
                    else filename
                    }


    input:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file(filtered_vcf), file(vcf_tbi) from filtered_data
    file genomefasta
    file genomefai
    file genomedict


    output:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file("${name}.biallelec.vcf.gz"), file("${name}.biallelec.vcf.gz.tbi") into BiallelecVCF

    script:

    """
    java -jar \$GATK_HOME/gatk-package-4.0.1.2-local.jar SelectVariants \\
    -R $genomefasta \\
    -V $filtered_vcf \\
    --restrict-alleles-to BIALLELIC \\
    -select-type SNP \\
    -O ${name}.biallelec.vcf
    bgzip -c ${name}.biallelec.vcf > ${name}.biallelec.vcf.gz 
    tabix -p vcf ${name}.biallelec.vcf.gz 
    """
}

 process allelespecificexpression {
    tag "$biallelec_vcf.baseName"
    publishDir "${params.outdir}/AlleleSpecificExpression", mode: 'copy'

    input:
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file(biallelec_vcf), file(vcf_index) from BiallelecVCF
    file genomefasta
    file genomefai
    file genomedict


    output:
    file "${name}.ASE.csv" into alleleSpecificExpression
    script:

    """
    java -jar \$GATK_HOME/gatk-package-4.0.1.2-local.jar ASEReadCounter \\
    -R $genomefasta \\
    -I $splitNCigar_bam \\
    -V $biallelec_vcf \\
    -O ${name}.ASE.csv
    """
}