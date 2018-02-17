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

# Install Fastqc
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o /opt/fastqc_v0.11.5.zip && \
    unzip /opt/fastqc_v0.11.5.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v0.11.5.zip

#install cutadapt for tirmgalore
RUN pip install cutadapt

# Install TrimGalore
RUN mkdir /opt/TrimGalore && \
    curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.2.zip -o /opt/TrimGalore/trim_galore_v0.4.2.zip && \
    unzip /opt/TrimGalore/trim_galore_v0.4.2.zip -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/trim_galore_v0.4.2.zip

# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# Install MultiQC
RUN pip install git+git://github.com/ewels/MultiQC.git

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

# Create root directories for common Swedish HPC systems
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw \
          /c3se /local /apps



