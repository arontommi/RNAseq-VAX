#!/usr/bin/env nextflow
/*

========================================================================================
                         RNAseq_VAX
========================================================================================
 RNAseq-VAX nalysis Pipeline. Started 2018-02-06.
 #### Homepage / Documentation
 RNAseq-VAX
 #### Authors
 Aron T. Skaftason arontommi <aron.skaftason@ki.se> - https://github.com/arontommi>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     RNAseq_VAX v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run RNAseq_VAX  -with-singularity rnaseq-vax.simg --project 'your_uppmax_project' \
    --fasta 'reference fasta used to align'
    
    Mandatory arguments:
      --fasta                      fasta used to align 
      --project                    your Uppmax project
      -with-singularity            Singularity container

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")


if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if (params.bamfolder ) {
    bam_md = Channel.fromPath(params.bamfolder+'*.bam')
}
else if ( !params.params.bamfolder ){
    exit 1, "No Bam specified! do some aligning first!"
}

if ( params.fasta ){
    if ( params.fasta.endsWith('.fa')) {
        genomefasta = file(params.fasta)
        genomefai = file(params.fasta + '.fai')
        genomedict = file(params.fasta - '.fa'+'.dict')
    }
    else if ( params.fasta.endsWith('.fasta')) {
        genomefasta = file(params.fasta)
        genomefai = file(params.fasta + '.fai')
        genomedict = file(params.fasta - '.fasta'+'.dict')
    }
}
/*
 * Readgroups added 
 */
process addReadGroups{ 
    tag "$bam_md.baseName"

    input:
    file bam_md

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
    picard AddOrReplaceReadGroups \\
        I= $bam_md \\
        O= ${name}.RG.bam \\
        RGLB=${params.addReadGroups.rglb} \\
        RGPL=${params.addReadGroups.rgpl} \\
        RGPU=${params.addReadGroups.rgpu} \\
        RGSM=${bam_md.baseName}

    samtools index ${bam_md.baseName}.RG.bam
    """
}
/*
 * SplitNCigarReads
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
    gatk SplitNCigarReads \\
    -R $genomefasta \\
    -I $rg_bam \\
    -O ${name}_split.bam 

    samtools index ${name}_split.bam
    """
}

/*
 * Haplotypecaller
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
    gatk HaplotypeCaller \\
    -R $genomefasta \\
    -I $splitNCigar_bam \\
    --dont-use-soft-clipped-bases \\
    --standard-min-confidence-threshold-for-calling ${params.haplotypeCaller.s_min_theshold} \\
    -O ${name}.vcf 
    bgzip -c ${name}.vcf > ${name}.vcf.gz
    tabix -p vcf ${name}.vcf.gz
    """
}

/*
 * varfiltering
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
    gatk VariantFiltration \\
    -R $genomefasta \\
    -V $vcf \\
    -window $params.varfiltering.window \\
    -cluster $params.varfiltering.cluster \\
    -filter-name FS \\
    -filter "FS > ${params.varfiltering.fs_filter}" \\
    -filter-name QD \\
    -filter "QD < ${params.varfiltering.qd_filter}"  \\
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
    set val(name), file("${name}.biallelec.vcf.gz"), file("${name}.biallelec.vcf.gz.tbi") into vep_annotating
    set val(name), file(splitNCigar_bam), file(splitNCigar_bam_bai), file("${name}.biallelec.vcf.gz"), file("${name}.biallelec.vcf.gz.tbi") into BiallelecVCF

    script:

    """
    gatk SelectVariants \\
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
    gatk ASEReadCounter \\
    -R $genomefasta \\
    -I $splitNCigar_bam \\
    -V $biallelec_vcf \\
    -O ${name}.ASE.csv
    """
}
