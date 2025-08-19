process FASTA_INDEX {
    tag "$meta.id"
    label "process_single"

    container "quay.io/jefffurlong/ref-utils"

    input:
    tuple val(meta), path(fasta)


    output:
    tuple val(meta), path("${fasta}"), emit: fasta
    tuple val(meta), path("${fasta.baseName}.dict"), emit: dict
    tuple val(meta), path("${fasta.baseName}.fasta.fai"), emit: fai

    script:
    """
    java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R=${fasta} O=${fasta.baseName}.dict
    samtools faidx ${fasta}

    """
}
