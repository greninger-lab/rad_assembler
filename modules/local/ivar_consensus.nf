process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/biocontainers/ivar:1.3.1--hecb563c_3"


    input:
        tuple val(meta), path(bam) 

    output:
        tuple val(meta), path("*_ivar_final_consensus.fasta"), emit: consensus
        tuple val(meta), path("*.bed"),                        emit: bed
    

    shell:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    
    # generate coverage bed file
    samtools depth -a -H ${bam} -o ${prefix}_final.bed

    # call consensus genome 
    samtools mpileup -d 5000 -A -Q 0 ${bam} | ivar consensus -p ${prefix} -n 'N' -m 5 -t 0 -i ${prefix}
    cp ${prefix}.fa ${prefix}_ivar_final_consensus.fasta

    """
}