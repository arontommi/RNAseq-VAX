FROM openjdk:8

LABEL authors=aron.skaftason@ki.se \
    description="Docker image containing all requirements for rnaseq-vax annotat pipeline"


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
        unzip \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o /opt/samtools-1.3.1.tar.bz2 && \
    tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.3.1 && \
    make && \
    make install && \
    rm /opt/samtools-1.3.1.tar.bz2

RUN conda install pandas click


ENV SNPEFF_VERSION=4_3r

# install snpSift
RUN mkdir /build && \
	cd /build && \
	wget --quiet -O snpEff_v${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_v${SNPEFF_VERSION}_core.zip \
  && unzip snpEff_v${SNPEFF_VERSION}_core.zip -d /opt/ \
  && rm snpEff_v${SNPEFF_VERSION}_core.zip