process READS_MAPPED {
    tag "$meta.id"
    label 'process_low'

    container "${ 'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), stdout,         emit: reads_mapped
    tuple val(meta), path("${meta.id}_empty.fasta"), emit: empty
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -F 4 -c ${bam}
    touch ${meta.id}_empty.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
