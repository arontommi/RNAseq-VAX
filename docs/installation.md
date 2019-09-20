# RNA-seq_VAX Installation

Currently this is only set up to run on Uppmax with singularity. 

first pull the singularity image like this :
```
singularity pull --name rnaseq_vax.sif docker://arontommi/rnaseq_vax:latest
```

Pull the project from github 
```
wget https://github.com/arontommi/RNAseq-VAX/archive/master.zip
mkdir -p ~/RNAseq-VAX
unzip master.zip -d ~/RNAseq-VAX/
cd ~/RNAseq-VAX/
```

And you are good go go! 