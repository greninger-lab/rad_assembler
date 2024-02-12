# base image
FROM python:3.7-bullseye

# install blast
RUN apt-get update \
    && apt-get install -y ncbi-blast+ python2.7 bwa bowtie2 samtools git autoconf wget build-essential \
    zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libsqlite3-dev libreadline-dev libffi-dev libbz2-dev liblzma-dev libncurses-dev nano default-jre

RUN wget https://www.python.org/ftp/python/3.7.16/Python-3.7.16.tar.xz
RUN tar -xf Python-3.7.16.tar.xz
WORKDIR /Python-3.7.16
RUN ./configure --enable-optimizations --enable-shared
RUN make -j 4
RUN make altinstall
RUN ldconfig /usr/local/share/python3.7

RUN /usr/local/bin/pip3.7 install biopython==1.81 numpy pandas

RUN apt-get install -y mummer libperl4-corelibs-perl bwa

WORKDIR /

RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.95.tar.gz && \
    tar -xzf BBMap_38.95.tar.gz && \
    rm BBMap_38.95.tar.gz

RUN git clone --depth 1 --branch v2.23 https://github.com/lh3/minimap2.git && (cd minimap2 && make)
WORKDIR /minimap2
RUN make

WORKDIR /

ADD scaffold_builder.py /scaffold_builder.py

ENV PATH=${PATH}:/usr/bin:/mugsy_x86-64-v1r2.2/MUMmer3.20:/minimap2:/bbmap

RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -f -p "$HOME"/miniconda
RUN rm Miniconda3-latest-Linux-x86_64.sh

RUN git clone -b v0.4.8 https://github.com/rrwick/Unicycler.git
WORKDIR Unicycler
RUN /usr/local/bin/python3.7 setup.py install

WORKDIR /

RUN wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
RUN tar -xzf SPAdes-3.13.0-Linux.tar.gz

ENV PATH=/SPAdes-3.13.0-Linux/bin:$PATH:/root/miniconda/bin

RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge  -c bioconda seqfu ivar picard biopython==1.81 bioconda::biopython.convert yuxiang::bam2fastq pilon

RUN update-alternatives --install /usr/bin/python python /usr/local/bin/python3.7 1

RUN /usr/local/bin/pip3.7 install pysam
# copy genomic id folders
ADD genome_identification /genome_identification
