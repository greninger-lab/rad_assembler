# rad_assembler
Nextflow pipeline for reference assisted de novo assembling and annotating whole genome shotgun sequenced viruses.

## Usage
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

Install [`Docker`](https://docs.docker.com/engine/installation/)

### Examples:<br>
    nextflow run greninger-lab/rad_assembler --input samplesheet.csv --outdir output_folder --genbank_ref NC_XXXXXXX.gb  -r main

#### To run it on AWS, add your nextflow config for aws after -c<br>
    nextflow run greninger-lab/rad_assembler --input samplesheet.csv --outdir output_folder --genbank_ref NC_XXXXXXXX.gb -profile docker -with-tower -c nextflow_aws.config -r main

#### Samplesheet example:<br>
assets/samplesheet.csv

##### You can create a samplesheet from an s3 folder (including single depth subdirectories) using the shell script generate_aws_samplesheet.sh in bin folder, like this:
generate_aws_samplesheet.sh s3://bucket-name /folder/path/to/run/folder/ samplesheet_name _L001_R1_001 _L001_R2_001

note:  replace _L001_R1_001 _L001_R2_001 with the suffixes of read1 and read2 if necessary
