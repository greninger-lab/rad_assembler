
FROM staphb/mummer:4.0.0-rgdv2-gnuplot

RUN apt-get update

RUN apt-get install -y samtools

RUN apt-get install -y bcftools

RUN apt-get install -y bowtie2

RUN apt-get install -y bwa

RUN apt-get install -y libboost-all-dev

RUN wget https://sourceforge.net/projects/gapfiller/files/v2.1.2/gapfiller-2.1.2.tar.gz && \
  tar -xvf gapfiller-2.1.2.tar.gz && \
  rm gapfiller-2.1.2.tar.gz

WORKDIR /data/gapfiller-2.1.2

RUN ./configure; make

ENV PATH=$PATH:/data/gapfiller-2.1.2/src

RUN cpan App::Cmd::Setup


WORKDIR /data

nucmer --maxmatch -b 4000 -c 30 -p W44547_pooneh NC_001806.fasta W44547_pooneh.scaffolds.fa
delta-filter -q W44547_pooneh.delta > W44547_pooneh.filter.out
show-tiling -c -R -p W44547_pooneh_result W44547_pooneh.filter.out > W44547_pooneh.tile.out

samtools mpileup –redo-BAQ –max-idepth 1000000 –max-depth 1000000 –min-MQ 20 –BCF –output-MQ –output-tags=INFO/AD,INFO/ADF,INFO/ADR,DP,SP
bcftools call –multiallelic-caller –variants-only
bcftools filter -e ’ALT=~.~ || QUAL<20 || AD[0]+AD[1]<100’



wget https://sourceforge.net/projects/samtools/files/samtools/0.1.17/samtools-0.1.17.tar.bz2
tar -xf samtools-0.1.17.tar.bz2
cd samtools-0.1.17
make