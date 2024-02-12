process GET_BEST_REFERENCE {
    tag "$meta.id"
    label "process_high"

    container "${ 'docker.io/jefffurlong/rad_utils' }"


    input:
    tuple val(meta), path(scaffolds), path(reads), path(mapped_ref_bam)
    val(find_ref_model)
    
    output:
    tuple val(meta), path("*.json"), emit: region_map
    tuple val(meta), path("region_*.fasta"), emit: regions
    tuple val(meta), path("*.gb"), emit: gb
    //tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.stats"), emit: stats
    //tuple val(meta), path("*.fasta")
    //tuple val(meta), path("*.bam")
    //tuple val(meta), path("*.sam")
    //tuple val(meta), path("*.fastq")

    script:
    """
    touch t.t
    /genome_identification/${find_ref_model}/get_top_hit.sh ${scaffolds}
    mapped_reads=`samtools view -F 4 -c ${mapped_ref_bam}`
    python3 /genome_identification/${find_ref_model}/extract_build_reference.py ${scaffolds} ${meta.id} ${reads[0]} ${reads[1]} \$mapped_reads ${task.cpus}

    """
}