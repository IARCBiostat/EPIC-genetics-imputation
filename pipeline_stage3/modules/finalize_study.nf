process FINALIZE_STUDY {
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${study_name}/stage3/final/${filename}"
    }

    input:
    tuple val(study_name), path(merged_pgen), path(merged_pvar), path(merged_psam), path(remove_list), path(sample_review_tsv)

    output:
    tuple val(study_name), path("${study_name}.pgen"), path("${study_name}.pvar"), path("${study_name}.psam")

    script:
    def merged_prefix = merged_pgen.baseName
    """
    if [ -s "${remove_list}" ]; then
      \$PLINK2_BIN \\
        --pfile ${merged_prefix} \\
        --remove ${remove_list} \\
        --make-pgen \\
        --out ${study_name}
    else
      \$PLINK2_BIN \\
        --pfile ${merged_prefix} \\
        --make-pgen \\
        --out ${study_name}
    fi
    """
}
