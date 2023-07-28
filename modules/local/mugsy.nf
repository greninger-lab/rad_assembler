process MUGSY {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'staphb/mugsy' }"

    input:
    tuple val(meta), path(scaffolds)
    path(fasta)

    output:
    tuple val(meta), path("*nonzero_ref.maf"), emit: aligned_maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # mugsy
    cp ${fasta} zzz_${fasta}
    mugsy --directory ./ --prefix ${meta.id}_aligned_scaffolds_ref zzz_${fasta} ${scaffolds}
    sed '/^a score=0/,\$d' ${meta.id}_aligned_scaffolds_ref.maf > ${meta.id}_aligned_scaffolds_nonzero_ref.maf

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mugsy: \$(echo 'mugsy 1r2.3')
    END_VERSIONS
    """
}
