# base image
FROM ubuntu:bionic

# install ubuntu dependencies
RUN apt-get update && apt-get install -y wget \
  build-essential \
  perl \
  perl-doc

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

