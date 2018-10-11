FROM nfcore/base
MAINTAINER Aron Skaftason <aron.skaftason@ki.se>
LABEL authors="aron.skaftason@ki.se" \
    description="Docker image containing all requirements for the RNAseq_VAX"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rnaseq-vax/bin:$PATH
