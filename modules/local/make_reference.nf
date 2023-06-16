process MAKE_REFERENCE {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'quay.io/michellejlin/tpallidum_wgs' }"

    input:
    tuple val(meta), path(sorted_bam)
    tuple val(meta), path(fasta)
    path(make_reference)

    output:
    tuple val(meta), path("*.fasta"), emit: new_ref
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # make_reference
    Rscript --vanilla make_reference.R ${sorted_bam} ${fasta} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        make_reference.R: \$(echo 'from cmv_make_reference.R Pavitra Roychoudhury Sep 2017')
    END_VERSIONS
    """
}
