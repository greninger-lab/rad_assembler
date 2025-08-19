process MUGSY {
    tag "$meta.id"
    label "process_high"

    container "quay.io/jefffurlong/mugsy_make_reference"

    input:
    tuple val(meta), path(scaffolds)
    path(fasta_ref)
    path(genbank_ref)
    path(regions)
    path(maf_convert)
    path(make_reference)
    path(region_map)
    //path(hybrid_ref)
    path(configure_reference)

    output:
    tuple val(meta), path("*.maf"), emit: aligned_maf
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*new_ref_consensus.fasta"), emit: new_ref
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id.replaceAll('-','_')
    def single_ref = regions ? 0 : 1

    """
    if [ ${single_ref} -eq 1 ]
    then
        # Use the entire reference
        ref_fasta="${fasta_ref}"
        ref_name=\${ref_fasta%.fasta*}
        mugsy --directory ./ --prefix ${prefix}_aligned_scaffolds_\${ref_name/./_} ${fasta_ref} ${scaffolds}
        sed '/^a score=0/,\$d' ${prefix}_aligned_scaffolds_\${ref_name/./_}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf
        /usr/bin/python2.7 ${maf_convert} -d sam ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam
        samtools view -bS -T ${fasta_ref} ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam | samtools sort > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam
        Rscript --vanilla make_reference.R ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam ${fasta_ref} 200
        mv ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}_consensus.fasta ${meta.id}_new_ref_consensus.fasta
    else
        # Loop through all regions of the reference and concatenate resulting new region references as defined in region_map
        array=( ${regions} )
        for ref in "\${array[@]}"
        do
            ref_name=\${ref%.fasta*}
            mugsy --directory ./ --prefix ${prefix}_aligned_scaffolds_\${ref_name/./_} \$ref ${scaffolds}
            sed '/^a score=0/,\$d' ${prefix}_aligned_scaffolds_\${ref_name/./_}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf
            minimumsize=500
            actualsize=\$(wc -c < "${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf")
            if [ \$actualsize -ge \$minimumsize ]; then
                /usr/bin/python2.7 ${maf_convert} -d sam ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam
                samtools view -bS -T \$ref ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam | samtools sort > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam
                Rscript --vanilla make_reference.R ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam \$ref 20
            else
                cp \$ref ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}_consensus.fasta
            fi
        done

        #cat *nonzero_region*consensus.fasta *nonzero_region*consensus.fasta > ${meta.id}_query.fasta
        #mv ${meta.id}_sb_Scaffold.fasta ${meta.id}_sb_new_ref_consensus.fasta

        /usr/bin/python3.6 configure_reference.py -r ${region_map} --build_reference -s ${meta.id}_aligned_scaffolds_nonzero_region ${genbank_ref}

        # mv ${meta.id}_aligned_scaffolds_nonzero_region_new_ref_consensus.fasta ${meta.id}_new_ref.fasta
        # mugsy --directory ./ --prefix ${prefix}_aligned_scaffolds_concat ${meta.id}_new_ref.fasta ${scaffolds}
        # sed '/^a score=0/,\$d' ${prefix}_aligned_scaffolds_concat.maf > ${prefix}_aligned_scaffolds_nonzero_concat.maf
        # /usr/bin/python2.7 ${maf_convert} -d sam ${prefix}_aligned_scaffolds_nonzero_concat.maf > ${prefix}_aligned_scaffolds_nonzero_concat.sam
        # samtools view -bS -T ${meta.id}_new_ref.fasta ${prefix}_aligned_scaffolds_nonzero_concat.sam | samtools sort > ${prefix}_aligned_scaffolds_nonzero_concat.bam
        # Rscript --vanilla make_reference.R ${prefix}_aligned_scaffolds_nonzero_concat.bam ${meta.id}_new_ref.fasta
        # mv ${prefix}_aligned_scaffolds_nonzero_concat_consensus.fasta ${meta.id}_sb_new_ref_consensus.fasta

    fi



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mugsy: \$(echo 'mugsy 1r2.2')
    END_VERSIONS
    """        
}
