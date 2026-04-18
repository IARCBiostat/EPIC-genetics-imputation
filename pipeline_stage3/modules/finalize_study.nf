process FINALIZE_STUDY {
    tag "${study_name}"
    cpus 2
    memory '16 GB'
    publishDir "${params.outdir}/${study_name}/stage3", mode: 'copy', pattern: "${study_name}.p*", overwrite: true

    input:
    tuple val(study_name), path(merged_pgen), path(merged_pvar), path(merged_psam), path(remove_list), path(sample_qc_tsv)

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
