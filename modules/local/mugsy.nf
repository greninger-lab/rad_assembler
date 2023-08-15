process MUGSY {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'docker.io/jefffurlong/mugsy:v1r2.2' }"

    input:
    tuple val(meta), path(scaffolds)
    path(fasta)

    output:
    tuple val(meta), path("*nonzero_ref.maf"), emit: aligned_maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id.replaceAll('-','_')
    """
    # mugsy
    mugsy --directory ./ --prefix ${prefix}_aligned_scaffolds_ref ${fasta} ${scaffolds}
    sed '/^a score=0/,\$d' ${prefix}_aligned_scaffolds_ref.maf > ${meta.id}_aligned_scaffolds_nonzero_ref.maf

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mugsy: \$(echo 'mugsy 1r2.3')
    END_VERSIONS
    """
}
