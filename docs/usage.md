# RNA-seq_VAX Usage

This pipeline is still in alpha stage. please use at your own cause


## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux`  inside a `.sh` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx6g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
    nextflow run RNAseq-VAX --genome 'The_the_ref_genome_used_to_align'  profile docker 
```

Yes it its that simple! ( if you have already run nf-core/rnaseq and have access to the internet) it looks for bam files in your result files and uses them for downstream processing.

another version of a run would be something like this:

```bash
    nextflow run RNAseq-VAX  -with-singularity rnaseq_vax.sif --project 'your_uppmax_project' --genome /sw/data/uppnex/reference/Homo_sapiens/hg19/program_files/GATK/concat.fasta
```


This will launch the pipeline with the `singularity` configuration profile ( how i have used it). generated from the docker image. A singularity image can be generated so: 

```bash
singularity build rnaseq_vax.sif docker://arontommi/rnaseq-vax:latest
```
you can also see that you can give a project id if you have to do so for your cluster.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull arontommi/RNAseq-VAX
```

