process BBMAP_DEDUPE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    dedupe.sh -Xmx\$maxmem threads=$task.cpus $raw out=${meta.id}.deduped.fq outd=${meta.id}.duplicates ac=f &> ${meta.id}.dedupe.log
    reformat.sh in=${meta.id}.deduped.fq out1=${meta.id}.deduped_1.fastq.gz out2=${meta.id}.deduped_2.fastq.gz
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
