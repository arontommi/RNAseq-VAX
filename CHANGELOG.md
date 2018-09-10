
## 0.1.0 - 2018-02-06
Initial release of RNA-seq_AVC, created with the NGI-NFcookiecutter template: https://github.com/ewels/NGI-NFcookiecutter.
## 0.1.1 - 2018-02-08
everything "stolen" and modified from [NGI-RNAseq](https://github.com/SciLifeLab/NGI-RNAseq)

dockerfile "stolen" and modified as well from same source

## 0.1.2 - 2018-03-20

first functional version, can take fastq files and mimic NGI-RNAseq or straight from bamfiles from final/markduplicates directory 
currently only works with star


## 0.1.3 - 2018-04-23

aim changed, now takes bamfiles from markDuplicates and runs the pipepline on that

## 0.1.4 - 2018-05-02

some code cleaning done

TODO:
	filter variants
	annotate variants 

## 0.1.5 - 2018-08-26

Renamed to RNAseq-VAX

annotation has been added. VEP has been added to main.nf and annotate.nf has also been created to create nice tab files. annotate.nf still needs to be cleaned up and fixed.

## 0.1.6 - 2018-09-05

strelka2pandas.py has been removed from bin and replaced by a similar app that works(vcf2tab project, will create a seperate github on that when it is worthy enough)

finally added a .gitignore file 

TODO: 
    fix annotation.nf
    create a seperated docker for annotation
    move vep to annotation ?
    set up so that i can run multiple dockers


## 0.1.7 - 2018-09-10

fixes made to vcf2tab
fixes made to annotation docker, made similar to og docker
fixes made to annotate.nf

filter added to Haplotypecaller2tab 
    needs at least depth of 10 alleles and 20 ref.


TODO:
    add resources to all processes of annoatate.nf into base.config
    see todo from 09-05
    fix igenomes to work
    


