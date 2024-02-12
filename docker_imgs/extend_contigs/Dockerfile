# base image
FROM python:bullseye

# install
RUN pip3 install biopython==1.81
RUN pip3 install numpy
RUN pip3 install pandas

# install blast
RUN apt-get update \
    && apt-get install -y bowtie2 samtools git

RUN git clone https://github.com/dzerbino/velvet
WORKDIR /velvet
RUN make

ADD extend_contigs.py /extend_contigs.py

ENV PATH=${PATH}:/velvet 


