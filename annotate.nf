#!/usr/bin/env nextflow
/*

========================================================================================
                         RNA-seq_AVC
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
params.vcf_dir = './results/VEP'
params.cosmic = false

cosmic = file(params.cosmic)

if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if (params.vcf_dir ) {
    vcfs = Channel.fromPath(params.vcf_dir+'/*.ann.vcf')

}

else if ( !params.params.vcf_dir ){
    exit 1, "No VCF specified! do some variant_calling first!"
}

process siftAddCosmic {
    tag {vcf}
    input:
       file vcf from vcfs
       file cosmic
    
    output:
        file("${vcf.baseName}.cosmic.ann.vcf") into siftAddCosmic

    script:
    """

    java -Xmx4g \
	  -jar /opt/snpEff/SnpSift.jar \
	  annotate $cosmic \
	  -info CNT \
	  $vcf \
	  > ${vcf.baseName}.cosmic.ann.vcf
    """
}

process finishVCF {
    tag {vcf}

	publishDir "${params.outdir}/csv_vcfs", mode: 'copy'
    
    input:
        file(vcf) from siftAddCosmic

    output:
        file("${vcf.baseName}.tab.csv") into finishedFile

    script:
    """

    python3 HaplotypeCaller2tab.py -i ${vcf} -o ${vcf.baseName}.tab.csv -s ${vcf.baseName}

    """ 

}
