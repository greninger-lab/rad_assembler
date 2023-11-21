process SRA_SCRUB_HUMAN {

    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'docker.io/jefffurlong/human_scrubber:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*scrubbed*.gz") , emit: reads


    script: 
    def fq1 = reads[0].getFileName().toString().replaceFirst(".gz", "")
    def fq2 = reads[1].getFileName().toString().replaceFirst(".gz", "")
    """

    gunzip ${reads[0]}
    gunzip ${reads[1]}

    seqfu ilv -1 ${fq1} -2 ${fq2} | /opt/scrubber/scripts/scrub.sh -x -i - -o - | seqfu dei -o ${meta.id}_scrubbed - 

    for file in *.fq; do
        mv -- "\$file" "\${file%.fq}.fastq"
        gzip "\${file%.fq}.fastq"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       ncbi/sra-human-scrubber:latest
    END_VERSIONS
    """
}
