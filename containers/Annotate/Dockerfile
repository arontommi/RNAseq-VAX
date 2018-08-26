FROM broadinstitute/gatk3:3.8-0

LABEL authors=aron.skaftason@ki.se \
    description="Docker image containing all requirements for rnaseq-vax pipeline"

# Export PATH
ENV \
GATK_HOME=/usr \
SNPEFF_VERSION=4_3r \
SNPEFF_HOME=/opt/snpEff \
PATH="/root/.pyenv/bin:${PATH}"

# Build container
ADD \
  build.sh /usr/bin/
RUN \
  build.sh