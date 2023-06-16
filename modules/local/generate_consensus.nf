process GENERATE_CONSENSUS {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'quay.io/michellejlin/tpallidum_wgs' }"

    input:
    tuple val(meta), path(sorted_bam)
    tuple val(meta), path(sorted_bam_bai)
    tuple val(meta), path(fasta)
    path(generate_consensus)

    output:
    tuple val(meta), path("${meta.id}_final_consensus.fasta"), emit: consensus
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # generate_consensus
    Rscript --vanilla generate_consensus.R ${meta.id} ${fasta} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_consensus.R: \$(echo 'from cmv_generate_consensus.R Pavitra Roychoudhury Sep 2017')
    END_VERSIONS
    """
}
