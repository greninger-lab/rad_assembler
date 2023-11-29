# rad_assembler
Nextflow pipeline for reference assisted de novo assembling whole genome shotgun sequenced viruses.

## Usage
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

Install [`Docker`](https://docs.docker.com/engine/installation/)


### Command line:
    nextflow run greninger-lab/rad_assembler \
        --input PATH_TO_SAMPLE_CSV \                      # required
        --outdir PATH_TO_OUTPUT_FOLDER \                  # required
        --bowtie2_host_index PATH_TO_HOST_BOWTIE2_INDEX \ # optional (path to bowtie2 index of host to use for filtering host DNA)
        --region_map PATH_TO_REGION_MAP_FILE \            # optional (path to region map json file)
        --genbank_ref ../NC_XXXXXXXX.gb \                 # required (genbank reference)
        --spades_flag \                                   # optional (default is meta, alternatives are careful or unicycler)
        -profile docker \                                 # required
        -with-tower \                                     # optional (use if you want to use Nextflow Tower)
        -c nextflow_aws.config \                          # optional (AWS account config info) 
        -r main                                           # required (use the github main branch)


#### Region map files:
An optional region map in json format for splitting the reference into sections so that a new reference for read mapping can be built more accurately.  This can be very helpful when there are large inverted repeat regions in the genome. 

Example region map file
-----------------------
    {
        "regions": {
            "TRL" : [0,9212],
            "UL" : [9213,117159],
            "IRL_IRS" : [117160,132604],
            "US" : [132605,145588],
            "TRS" : [145589,152221]
        },
        "region_order": [
            {"region": "TRL", "reverse": false}, 
            {"region": "UL", "reverse": false}, 
            {"region": "IRL_IRS", "reverse": false}, 
            {"region": "US", "reverse": false},
            {"region": "TRS", "reverse": false}
        ]
    }
-----------------------

"regions" defines a region name and the reference coordinates for extracting the region.  De novo assembled scaffolds are mapped to the region and then a new reference is generated for the region. 
"region_order" sets the order for concatenating the new region references together to make a complete genome reference.  "reverse" will reverse compliment the region reference before concatenating if set to true.

#### Example region map files:
    rad_assembler/region_maps/HSV1-NC001806.json
    rad_assembler/region_maps/HSV2-NC001798.json


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
