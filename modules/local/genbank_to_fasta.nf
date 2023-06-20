process GENBANK_TO_FASTA {

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'biocontainers/biopython:1.75' }"

    input:
    path(genbank_file)

    output:
    path("${genbank_file.baseName}.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "from Bio import SeqIO; SeqIO.convert('${genbank_file}', 'genbank', '${genbank_file}'.replace('.gb', '.fasta'), 'fasta');"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython 1.75
    END_VERSIONS
    """
}


/*
process GENBANK_TO_FASTA {

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'biocontainers/biopython:1.75' }"

    input:
    path(genbank_file)

    output:
    path "*.fasta"       , emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    """
    python3 -c "from Bio import SeqIO; SeqIO.convert(${genbank_file}, genbank, ${genbank_file%.*}.fasta, fasta);"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
*/