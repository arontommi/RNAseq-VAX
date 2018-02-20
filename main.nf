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

    nextflow run rna-seq_AVC --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      --gtf                         Path to GTF file
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}



params.reads = false
params.deduped_bam = false 
params.singleEnd = false
params.outdir = './results'
params.genome = false
params.project = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.download_fasta = false
params.download_gtf = false


params.rglb = '1'
params.rgpl = 'lib1'
params.rgpu = 'illumina'


wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.saveTrimmed = false


// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
//forward_stranded = params.forward_stranded
//reverse_stranded = params.reverse_stranded
//unstranded = params.unstranded

// hisat2 has been removed (STAR masterrace (for now))

params.aligner = 'star'

// Validate inputs
if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeBED12;
              gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
}
else if ( !params.download_gtf ){
    exit 1, "No GTF annotation specified!"
}

if ( params.fasta ){
    if ( params.fasta.endsWith('.fa')) {
        fasta = file(params.fasta)
        fai = file(params.fasta + '.fai')
        dict = file(params.fasta - '.fa'+'.dict')
    }
    else if ( params.fasta.endsWith('.fasta')) {
        fasta = file(params.fasta)
        fai = file(params.fasta + '.fai')
        dict = file(params.fasta - '.fasta'+'.dict')
    }
}

if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if ( params.reads) { 
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { 
            exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_files_fastqc; read_files_trimming }
}
else if (params.deduped_bam) {
    Channel
        .fromFilePairs('${params.outdir}/markDuplicates/*.{bam,bam.bai}')
        .println{}
}