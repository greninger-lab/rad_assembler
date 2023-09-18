process GENBANK_TO_FASTA {

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'biocontainers/biopython:1.75' }"

    input:
    path(genbank_file)
    path(region_map)
    path(configure_reference)

    output:
    path("${genbank_file.baseName}.fasta"), emit: fasta
    path("region_*.fasta"), emit: regions, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def regionmap = region_map ? "-r ${region_map}" : ""
    """
    python3 configure_reference.py ${regionmap} ${genbank_file}
 
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