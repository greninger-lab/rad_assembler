process SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/':
        'staphb/samtools:1.17' }"

    input:
    tuple val(meta), path(rlog), path(qlog), path(align_new_ref_bam), path(align_new_ref_bai), path(consensus)

    output:
    path("*.tsv"), emit: summary_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # raw reads and trimmed reads  
    raw_reads=`grep "Input:" ${rlog} | awk '{print \$2}'`
    trimmed_reads=`grep "Result:" ${qlog} | awk '{print \$2}'` 
    
    pct_reads_trimmed=\$(python3 -c "print (float('\$trimmed_reads') / float('\$raw_reads') * 100)")

    # mapped reads
    mapped_reads=`samtools view -F 4 -c ${align_new_ref_bam}`
    pct_reads_mapped=\$(python3 -c "print (float('\$mapped_reads') / float('\$raw_reads') * 100)")
    # whole genome coverage
    pct_genome_covered=`samtools coverage ${align_new_ref_bam} | awk 'NR>1' | cut -f6`
    mean_genome_coverage=`samtools coverage ${align_new_ref_bam} | awk 'NR>1' | cut -f7`

    # consensus genome
    consensus_length=0
    num_ns_consensus=0
    num_as_consensus=0
    num_cs_consensus=0
    num_gs_consensus=0
    num_ts_consensus=0
    num_non_ns_ambiguous=0
    
    filesize=$(wc -c < ${consensus})
    if ((filesize > 500)); then
        consensus_length=`awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${consensus} | awk 'FNR==2{print val,\$1}'`
        num_ns_consensus=`grep -v "^>" ${consensus} | tr -c -d N | wc -c`
        num_as_consensus=`grep -v "^>" ${consensus} | tr -c -d A | wc -c`
        num_cs_consensus=`grep -v "^>" ${consensus} | tr -c -d C | wc -c`
        num_gs_consensus=`grep -v "^>" ${consensus} | tr -c -d G | wc -c`
        num_ts_consensus=`grep -v "^>" ${consensus} | tr -c -d T | wc -c`
        num_non_ns_ambiguous=\$(python3 -c "print(int(\$consensus_length)-int(\$num_as_consensus)-int(\$num_cs_consensus)-int(\$num_gs_consensus)-int(\$num_ts_consensus)-int(\$num_ns_consensus))")
    fi

    echo "sample_name\traw_reads\ttrimmed_reads\tpct_reads_trimmed\tmapped_reads\tpct_reads_mapped\tpct_genome_covered\tmean_genome_coverage\tconsensus_length\tnum_ns\tnum_ambiguous" > ${prefix}.summary.tsv
    echo "${prefix}\t\${raw_reads}\t\${trimmed_reads}\t\${pct_reads_trimmed}\t\${mapped_reads}\t\${pct_reads_mapped}\t\${pct_genome_covered}\t\${mean_genome_coverage}\t\${consensus_length}\t\${num_ns_consensus}\t\${num_non_ns_ambiguous}" >> ${prefix}.summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summary: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
