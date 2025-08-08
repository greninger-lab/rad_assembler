# base image
FROM quay.io/fedora/python-312:312

# install biocode
RUN pip3 install biopython==1.81
RUN pip3 install biocode
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install biopython-convert

USER root

RUN dnf update -y && \
    dnf install -y samtools java-17-openjdk && \
    dnf clean all && \
    rm -rf /var/cache/dnf

RUN curl -L -o /usr/local/bin/picard.jar \
    https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

