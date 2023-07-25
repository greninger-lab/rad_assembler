process BWA_MEM_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "dukegcb/bwa-samtools"


    input: 
        tuple val(meta), path(reads)  //, file("${base}_summary.csv")
        tuple val(meta2), path(reference_fasta)
    output:
        tuple val(meta), path("*_new_ref.bam"), emit: new_ref_bam
        tuple val (meta), path("*.txt"),        emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"    
    """
    #!/bin/bash

    /usr/local/bin/bwa index ${reference_fasta}
    /usr/local/bin/bwa mem -t ${task.cpus} ${reference_fasta} $raw | samtools view -@ ${task.cpus} -b -F 4 - > ${prefix}_new_ref.bam
    reads_mapped=\$(samtools view -c ${base}_new_ref.bam)

    #cp ${prefix}_summary.csv ${prefix}_summary2.csv
    touch ${prefix}_new_ref_reads_mapped.txt
    printf "\$reads_mapped" >> ${prefix}_new_ref_reads_mapped.txt

    """
}