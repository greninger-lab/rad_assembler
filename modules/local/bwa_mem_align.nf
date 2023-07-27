process BWA_MEM_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/dukegcb/bwa-samtools"


    input: 
        tuple val(meta), path(reads)  //, file("${base}_summary.csv")
        tuple val(meta2), path(reference_fasta)
    output:
        tuple val(meta), path("*_sorted.bam"), emit: new_ref_bam
        tuple val (meta), path("*.txt"),        emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"    
    """
    #!/bin/bash

    /usr/local/bin/bwa index ${reference_fasta}
    /usr/local/bin/bwa mem -t ${task.cpus} ${reference_fasta} $raw | samtools view -S -@ ${task.cpus} -b -F 4 | samtools sort -@ ${task.cpus} - > ${prefix}_sorted.bam
    reads_mapped=\$(samtools view -c ${prefix}_sorted.bam)

    #cp ${prefix}_summary.csv ${prefix}_summary2.csv
    touch ${prefix}_reads_mapped.txt
    printf "\$reads_mapped" >> ${prefix}_reads_mapped.txt

    """
}