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
        Rscript --vanilla make_reference.R ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam ${fasta_ref}
        mv ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}_consensus.fasta ${meta.id}_new_ref_consensus.fasta
    else
        # Loop through all regions of the reference and concatenate resulting new region references as defined in region_map
        array=( ${regions} )
        for ref in "\${array[@]}"
        do
            ref_name=\${ref%.fasta*}
            mugsy --directory ./ --prefix ${prefix}_aligned_scaffolds_\${ref_name/./_} \$ref ${scaffolds}
            sed '/^a score=0/,\$d' ${prefix}_aligned_scaffolds_\${ref_name/./_}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf
            /usr/bin/python2.7 ${maf_convert} -d sam ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.maf > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam
            samtools view -bS -T \$ref ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.sam | samtools sort > ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam
            Rscript --vanilla make_reference.R ${meta.id}_aligned_scaffolds_nonzero_\${ref_name}.bam \$ref 
        done

        python3 configure_reference.py -r ${region_map} --build_reference -s ${meta.id}_aligned_scaffolds_nonzero_region ${genbank_ref}
    fi



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mugsy: \$(echo 'mugsy 1r2.2')
    END_VERSIONS
    """        
}
