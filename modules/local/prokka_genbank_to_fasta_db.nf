process PROKKA_GENBANK_TO_FASTA_DB {

    label 'process_low'

    conda "bioconda::prokka=1.14.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka%3A1.14.6--pl5321hdfd78af_4' :
        'biocontainers/prokka:1.14.6--pl5321hdfd78af_4' }"

    input:
    path(genbank_file)

    output:
    path("${genbank_file.baseName}.faa"), emit: faa
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prokka-genbank_to_fasta_db ${genbank_file} > ${genbank_file.baseName}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
