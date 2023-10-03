process ABYSS {
    tag "$scaffolds"
    label 'process_high'

    //conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abyss:2.1.0' :
        'docker.io/pegi3s/abyss:2.1.0' }"


    input:
    tuple val(meta), path(scaffolds), path(reads)

    output:
    tuple val(meta), path("*_sealer_scaffold.fa"), emit: sealed_scaffolds
    tuple val(meta), path("*.txt"), emit: log
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    abyss-sealer -P10 -k90 -k80 -k70 -k60 -k50 -k40 -k30 -b500M -o ${meta.id}_sealer -S ${scaffolds} ${reads[0]} ${reads[1]}    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abyss: abyss:v2.1.5-7
    END_VERSIONS
    """

}
