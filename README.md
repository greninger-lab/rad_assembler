# rad_assembler
Nextflow pipeline for reference assisted de novo assembling whole genome shotgun sequenced viruses.

## Usage
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

Install [`Docker`](https://docs.docker.com/engine/installation/)


### Example command line usage:
    nextflow run greninger-lab/rad_assembler -r find_reference -latest \
              --input sample_fastqs.csv \
              --outdir ./out/ \
              --genbank_ref NC_001798.gb \
              --find_reference hsv2 \
              --spades_flag unicycler \
              --kraken_host_db s3://fh-pi-jerome-k-eco/greninger-lab/greninger-lab-file-share/refs/Kraken2_human/k2_human/ \
              -profile docker \
              -c ~/nextflow_aws.config \

## Command line options
| option | description | 
|--------|-------------|
| `--input  /path/to/sample_fastqs.csv` | (required) path to a csv sample,fastq_1,fastq_2 input file |
| `--outdir /path/to/output`                | (required) output directory |
| `--kraken_host_db /path/to/kraken2_human_db`  | (required) path to Kraken2 human database |
| `--genbank_ref <file>`        | (required) path to a GenBank (.gb) format reference file| 
| `--spades_flag <key>`        | (optional) default key is "meta", alternatives are "careful" or "unicycler" |
| `--find_reference <key>`     | (optional) available keys are "hsv1", "hsv2" or "cmv" |
| `-profile docker`                         | (required) |
| `-c /path/to/your/custom.config`          | (optional) used specify a custom configuration file (see [Nextflow docs](https://www.nextflow.io/docs/latest/config.html) |


#### Sample csv example:
    rad_assembler/assets/example.csv

#### Sample csv format:
---------
    sample,fastq_1,fastq_2
    F79217_S26,sample_L001_R1_001.fastq.gz,sample_L001_R2_001.fastq.gz
---------

You can create a sample csv file from an s3 folder (including single depth subdirectories) using the shell script generate_aws_sample_csv.sh in bin folder, like this:

    generate_aws_sample_csv.sh s3://bucket-name /folder/path/to/run/folder/ csv_name _L001_R1_001 _L001_R2_001
note:  replace _L001_R1_001 _L001_R2_001 with the suffixes of read1 and read2 if necessary
