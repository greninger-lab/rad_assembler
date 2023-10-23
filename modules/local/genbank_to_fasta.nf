process GENBANK_TO_FASTA {

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'docker.io/jefffurlong/biopython_biocode' }"

    input:
    path(genbank_file)
    path(region_map)
    path(configure_reference)
    path(convert_genbank_to_gff)

    output:
    path("${genbank_file.baseName}.fasta"), emit: fasta
    path("region_*.fasta"), emit: regions, optional: true
    path("${genbank_file.baseName}.gff"), emit: genes_gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def regionmap = region_map ? "-r ${region_map}" : ""
    """
    python3 configure_reference.py ${regionmap} ${genbank_file}
    python3 convert_genbank_to_gff3.py -i ${genbank_file} --no_fasta -o ${genbank_file.baseName}.gff
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert_genbank_to_gff3.py
    END_VERSIONS
    """
}
