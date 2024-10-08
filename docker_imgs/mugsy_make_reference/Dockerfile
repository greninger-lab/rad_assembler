# base image
FROM quay.io/bedrock/ubuntu:bionic-20230530

# install ubuntu dependencies
RUN apt-get update && apt-get install -y wget \
   build-essential \
   perl \
   perl-doc \
   python2.7

# install mugsy
RUN wget https://sourceforge.net/projects/mugsy/files/mugsy_x86-64-v1r2.2.tgz && \
  tar -xvf mugsy_x86-64-v1r2.2.tgz && \
  rm mugsy_x86-64-v1r2.2.tgz

# install perl dependency
RUN perl -MCPAN -e install Pod::Usage

# add path to mugsyenv.sh
RUN cd mugsy_x86-64-v1r2.2/ && \
  bash -c "sed -i 's/\.\/mugsy_x86-64-v1r2.3/\/mugsy_x86-64-v1r2.2/g' mugsyenv.sh"

RUN mkdir /data
WORKDIR /data

# export paths
ENV MUGSY_INSTALL="/mugsy_x86-64-v1r2.2"
ENV PATH="/mugsy_x86-64-v1r2.2:/mugsy_x86-64-v1r2.2/mapping:${PATH}"
ENV PERL5LIB="/mugsy_x86-64-v1r2.2/perllibs"

#####################
# add essential tools
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		software-properties-common \
        ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
        && add-apt-repository --enable-source --yes "ppa:marutter/rrutter3.5" \
	&& add-apt-repository --enable-source --yes "ppa:marutter/c2d4u3.5" 

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## This was not needed before but we need it now
ENV DEBIAN_FRONTEND noninteractive

# Now install R and littler, and create a link for littler in /usr/local/bin
# Default CRAN repo is now set by R itself, and littler knows about it too
# r-cran-docopt is not currently in c2d4u so we install from source
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
         littler \
 		 r-base \
 		 r-base-dev \
 		 r-recommended \
  	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
 	&& install.r docopt \
 	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
 	&& rm -rf /var/lib/apt/lists/*

CMD ["bash"]

## Install rJava
RUN apt-get -y update && apt-get install -y \
   default-jdk  r-cran-rjava  r-cran-nloptr libssh2-1-dev

RUN yes | apt-get install libv8-dev
RUN yes | apt-get install libxml2-dev
RUN yes | apt-get install libudunits2-dev
RUN yes | apt install libgdal-dev
RUN yes | apt-get install libmagick++-dev

## Install extra R packages using requirements.R
RUN R -e "install.packages('BiocManager'); library('BiocManager')"
RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('GenomicAlignments')"
RUN R -e "install.packages('Biostrings',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('RCurl')"
RUN R -e "BiocManager::install('doParallel')"
RUN R -e "BiocManager::install('ShortRead')"

SHELL ["/bin/bash", "-c"]
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.0-3-Linux-x86_64.sh && \
	bash Miniconda3-py38_23.5.0-3-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-py38_23.5.0-3-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc"
RUN conda install -c bioconda bowtie2 bedtools bwa mafft bcftools tabix
RUN conda clean -afy

WORKDIR /usr/src
#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

RUN apt-get install -y python3-pip
RUN pip3 install biopython

ENV PATH=${PATH}:/usr/src/samtools-1.9:/usr/bin
