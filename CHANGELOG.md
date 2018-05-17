
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
	- add annotation []
		- dbsnp within haplotypecaller [x]
		- cosmic []
		- add gene annotation for vcf files []
	- fix if statements for different parameters so that error messages are presented when files are not found or given []
		- for if params.bamfolder []
 		- for if params.genome []
 	- fix channel input from an if statement []
 	- clean up steps for each process []
 	- fix containers so that each process has its own container []
 	- fix dbsnp and gene stuff []
 	- make sure annotate_hpc process works
 	- add tests []
 	
 	- finish todo list :P []








