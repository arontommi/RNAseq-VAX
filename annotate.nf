#!/usr/bin/env nextflow
/*

========================================================================================
                         ANNOTATE_VAX
========================================================================================
 RNA-seq_AVC Analysis Pipeline. Started 2018-02-06.
 #### Homepage / Documentation
 RNAseq-VAX
 #### Authors
 Aron T. Skaftason arontommi <aron.skaftason@ki.se> - https://github.com/arontommi>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     RRNAseq_VAX.annotate v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run RNAseq_VAX/annotate.nf -with-singularity annotate.img --project 'your_uppmax_project' --cosmic /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/GRCh38_b37_cosmic_v74.noCHR.vcf

    Mandatory arguments:
      --cosmic                      what cosmic file to annotate with
      --project                     your Uppmax project
      -with-singularity             Singularity container


    """.stripIndent()
}

params.outdir = './results'
params.vcf_dir = './results/VariantFiltration'
params.cosmic = false
params.genome = 'GRCh37'
params.vep_dir = false
params.vep_version = 95 


cosmic = file(params.cosmic)

if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if (params.vcf_dir ) {
    vcfs = Channel.fromPath(params.vcf_dir+'/*.vcf.gz')

}

else if ( !params.params.vcf_dir ){
    exit 1, "No VCF specified! do some variant_calling first!"
}

process RunVEP {
    tag "$vcffile.baseName"
    publishDir "${params.outdir}/VEP", mode: 'copy'
    println "container : $workflow.container"
    println "profile : $workflow.profile"

    input:
    file(vcffile) from vcfs
    val genome from params.genome
    val vep_version from params.vep_version
    output:
    set val(name), file("${name}.VEP.summary.html"), file("${name}.VEP.ann.vcf") into vep_out

    script:
    name = "${vcffile.baseName}"

    """
    vep --dir /.vep \
        -i $vcffile \
        -o ${name}.VEP.ann.vcf \
        --assembly $genome \
        --cache \
        --cache_version $vep_version \
        --database \
        --everything \
        --filter_common \
        --format vcf \
        --offline \
        --per_gene \
        --stats_file ${name}.VEP.summary.html \
        --total_length \
        --vcf
    """
}

process siftAddCosmic {
    tag {vcf}
    input:
        set val(name), file(html), file(vcf) from vep_out
        file cosmic
    
    output:
        file("${vcf.baseName}.cosmic.ann.vcf") into siftAddCosmic

    script:
    """

    SnpSift annotate $cosmic \
	  -info CNT \
	  $vcf \
	  > ${vcf.baseName}.cosmic.ann.vcf
    """
}
