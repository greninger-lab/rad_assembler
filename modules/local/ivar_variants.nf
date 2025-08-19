process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/jefffurlong/ivar:1.4.4"


    input:
        tuple val(meta), path(bam)
        path(ref)
        path(gff)


    output:
        tuple val(meta), path("*.tsv"), emit: variants
    

    shell:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash

    # call variants
    samtools \\
        mpileup \\
        --ignore-overlaps --count-orphans --no-BAQ --max-depth 0 --min-BQ 0 \\
        --reference ${ref} \\
        ${bam} \\
        | ivar \\
            variants \\
            -q ${params.ivar_variants_q} -t ${params.ivar_variants_t} -m ${params.ivar_variants_m} -G \\
            -g ${gff} \\
            -r ${ref} \\
            -p ${prefix}

    """
}
