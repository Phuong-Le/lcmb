FROM ubuntu:24.04
LABEL maintainer="Phuong Le" 
LABEL maintainer.email="al35@sanger.ac.uk"


# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# install python and conda 
RUN apt-get update \
    && apt-get upgrade -y \
	&& apt-get install -yq \
        build-essential \
        wget \
        default-jdk \
	&& wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh \
	&& /bin/bash ./miniconda.sh -b -p /opt/conda \
	&& rm -f ./miniconda.sh \
	&& /opt/conda/bin/conda init bash \
    && apt-get -y install git-all \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
ENV PATH=/opt/conda/bin:$PATH

# installing stan 

# installing R packages
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && conda install -y conda-forge::r-base=4.3.3 \
        conda-forge::r-common \
        conda-forge::r-devtools \ 
        conda-forge::r-codetools \
        conda-forge::r-biocmanager \
        conda-forge::r-ggplot2 \
        conda-forge::r-data.table \
        conda-forge::r-ggpubr \
        conda-forge::r-logr \
        conda-forge::r-bh \
        conda-forge::r-rstan \
        conda-forge::r-loo \
        conda-forge::r-scales \
        conda-forge::r-optparse 
RUN R -e "devtools::install_github('Phuong-Le/distributionsR')"

# testing that all R packages are installed correctly
RUN R -e "library(ggplot2); \
    library(data.table); \
    library(ggpubr); \
    library(logr); \
    library(rstan); \
    library(loo); \
    library(scales); \
    library(optparse); \
    library(distributionsR); \
    print('all packages are ready to be used')"
