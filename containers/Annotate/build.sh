#!/bin/sh
set -eu pipefail

apt-get update && apt-get install -y --no-install-recommends \
  python \
  python3 \
  git \
  python3-pip \
  python-pip \
  python-dev \
  python3-dev \
  libatlas-base-dev \
  gfortran \
  python3-pandas \
  python-numpy

# Install tools
mkdir /build

# install PyVCF - required by seqtool
cd /build
git clone https://github.com/jamescasbon/PyVCF.git
cd PyVCF
python setup.py install

# install seqtool
cd /build
git clone https://github.com/seandavi/SDST.git
cd SDST
python setup.py install

# install Python 3.6.3 and some fixes. NB - Python 3.6.3 globaled from here on and down
cd /build
apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
xz-utils tk-dev libffi-dev
curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv install 3.6.3
pyenv global 3.6.3
pip install click
pip install pandas

# install snpSift
cd /build
wget --quiet -O snpEff_v${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_v${SNPEFF_VERSION}_core.zip \
  && unzip snpEff_v${SNPEFF_VERSION}_core.zip -d /opt/ \
  && rm snpEff_v${SNPEFF_VERSION}_core.zip