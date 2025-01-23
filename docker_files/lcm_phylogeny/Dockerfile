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
    && apt-get install -y python3-pip python3-dev \
	&& wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh \
	&& /bin/bash ./miniconda.sh -b -p /opt/conda \
	&& rm -f ./miniconda.sh \
	&& /opt/conda/bin/conda init bash \
    && apt-get -y install git-all \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
ENV PATH=/opt/conda/bin:$PATH

# install pandas 
RUN pip install pandas==1.5.3 \
    && pip show pandas

# downloading mpboot 
RUN cd /opt \
    && wget http://www.iqtree.org/mpboot/mpboot-avx-1.1.0-Linux.tar.gz \
    && tar xvzf mpboot-avx-1.1.0-Linux.tar.gz \
    && cd mpboot-avx-1.1.0-Linux/bin \
    && mv mpboot-avx mpboot \
    && mv mpboot /usr/local/bin \
    && mpboot --help 

# installing R
RUN conda install -y -c conda-forge \
        r-base=4.3.3 \
        r-devtools \ 
        r-biocmanager \
        r-ape \
        r-ggplot2 
RUN conda install -y -c bioconda \
        bioconductor-ggtree 
RUN R -e "devtools::install_git('https://github.com/NickWilliamsSanger/treemut')"
# testing that all R packages are installed correctly
RUN R -e "library(ape); \
    library(ggplot2); \
    library(ggtree); \
    library(treemut); \
    print('all packages are ready to be used')"
