FROM openjdk:8

LABEL authors="{{ cookiecutter.author_email }}" \
    description="Docker image containing all requirements for {{ cookiecutter.pipeline_name }} pipeline"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libgsl2 \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        make \
        python-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py


# Install GATK
RUN curl -fsSL https://github.com/broadinstitute/gatk/releases/download/4.0.1.2/gatk-4.0.1.2.zip -o /opt/gatk-4.0.1.2.zip && \
    unzip /opt/gatk-4.0.1.2.zip -d /opt/ && \
    rm /opt/gatk-4.0.1.2.zip
ENV GATK_HOME /opt/gatk-4.0.1.2

# Install PicardTools
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.0.1/picard-tools-2.0.1.zip -o /opt/picard-tools-2.0.1.zip && \
    unzip /opt/picard-tools-2.0.1.zip -d /opt/ && \
    rm /opt/picard-tools-2.0.1.zip
ENV PICARD_HOME /opt/picard-tools-2.0.1

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o /opt/samtools-1.3.1.tar.bz2 && \
    tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.3.1 && \
    make && \
    make install && \
    rm /opt/samtools-1.3.1.tar.bz2

# Install bcftools
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/1.7/bcftools-1.7.tar.bz2 -o /opt/bcftools-1.7.tar.bz2 && \
    tar xvjf /opt/bcftools-1.7.tar.bz2 -C /opt/ && \
    cd /opt/bcftools-1.7 && \
    make && \
    make install && \
    rm /opt/bcftools-1.7.tar.bz2

# Install htslib
RUN curl -fsSL https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 -o /opt/htslib-1.7.tar.bz2 && \
    tar xvjf /opt/htslib-1.7.tar.bz2 -C /opt/ && \
    cd /opt/htslib-1.7 && \
    make && \
    make install && \
    rm /opt/htslib-1.7.tar.bz2


# Create root directories for common Swedish HPC systems
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw \
          /c3se /local /apps

# Setup ENV variables
ENV \
  SNPEFF_VERSION=4_3r \
  SNPEFF_HOME=/opt/snpEff

# Install snpEff
RUN \
  wget --quiet -O snpEff_v${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_v${SNPEFF_VERSION}_core.zip \
  && unzip snpEff_v${SNPEFF_VERSION}_core.zip -d /opt/ \
  && rm snpEff_v${SNPEFF_VERSION}_core.zip

ENV GENOME="GRCh37.75"

# Download Genome
RUN java -jar $SNPEFF_HOME/snpEff.jar download $GENOME

