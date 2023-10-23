
process FORMAT_VARIANTS {

    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'docker.io/jefffurlong/biopython_biocode' }"

    input:
    tuple val(meta), path(ivar_variants)
    path(gff)
    path(edit_ivar_variants)

    output:
    tuple val(meta), path("${meta.id}.formatted_variants.tsv"), emit: formatted_variants
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    """
    python3 ${edit_ivar_variants} ${ivar_variants} ${gff} ${meta.id}.formatted_variants
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edit_ivar_variants.py
    END_VERSIONS
    """
}
