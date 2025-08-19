
process FORMAT_VARIANTS {
    tag "$meta.id"
    label "process_single"

    container "quay.io/jefffurlong/ref-utils"

    input:
    tuple val(meta), path(ivar_variants)
    path(ref)
    path(gff)
    path(edit_ivar_variants)

    output:
    tuple val(meta), path("${meta.id}.formatted_variants.tsv"), emit: formatted_variants
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    """
    python3 ${edit_ivar_variants} ${ivar_variants} ${gff} ${ref} ${meta.id}.formatted_variants
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edit_ivar_variants.py
    END_VERSIONS
    """
}
