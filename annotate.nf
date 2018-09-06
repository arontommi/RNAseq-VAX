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
     RNA-seq_AVC.annotate v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run RNAseq_AVC/annotate.nf -with-singularity annotate.img --project 'your_uppmax_project' --cosmic '/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/GRCh38_b37_cosmic_v74.noCHR.vcf'

    Mandatory arguments:
      --cosmic                      what cosmic file to annotate with
      --project                     your Uppmax project
      -with-singularity             Singularity container


    """.stripIndent()
}

params.outdir = './results'
params.vcf_dir = './results/VEP'
params.cosmic = false


if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-modules' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

if (params.vcf_dir ) {
    Channel
        .fromPath(params.vcf_dir+'*.ann.vcf.gz')
        .set{vcfs}
}

else if ( !params.params.vcf_dir ){
    exit 1, "No VCF specified! do some variant_calling first!"
}

process siftAddCosmic {
    tag {vcf}
    input:
       file(vcf) from vcfs
       set file(cosmic) from params.cosmic
    
    output:
        file("${vcf.baseName}.cosmic.ann.vcf") into filteredcosmicvcf

    script:
    """

    grep -E '#|PASS' ${vcf} > ${vcf.baseName}.pass.vcf

    java -Xmx4g \
	  -jar \$SNPEFF_HOME/SnpSift.jar \
	  annotate \
	  -info CNT \
    ${cosmic} \
	  ${vcf.baseName}.pass.vcf \
	  > ${vcf.baseName}.cosmic.ann.vcf
    """
}

process finishVCF {
    tag {vcf}

	publishDir "${params.outdir}/AnnotatedFilteredVcfs", mode: 'copy'
    
    input:
        file(vcf) from siftAddCosmic

    output:
        file("${vcf.baseName}.tab.csv") into finishedFile

    script:
    """

    python HaplotypeCaller2tab.py -i ${vcf} -o ${vcf.baseName}.tab.csv -s ${vcf.baseName}

    """ 

}
