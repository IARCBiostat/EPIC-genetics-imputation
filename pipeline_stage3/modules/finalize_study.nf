process FINALIZE_CHROM {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${study_name}/stage3/final/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(pgen), path(pvar), path(psam), path(remove_list), path(sex_update), path(pheno_update)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.pgen"), path("${study_name}_chr${chr}.pvar"), path("${study_name}_chr${chr}.psam")

    script:
    def pfile_prefix = pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 120
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    if [ -s "${remove_list}" ]; then
      \$PLINK2_BIN \\
        --pfile ${pfile_prefix} \\
        --update-sex ${sex_update} \\
        --pheno ${pheno_update} \\
        --remove ${remove_list} \\
        --threads ${threads} \\
        --make-pgen \\
        --out "${study_name}_chr${chr}"
    else
      \$PLINK2_BIN \\
        --pfile ${pfile_prefix} \\
        --update-sex ${sex_update} \\
        --pheno ${pheno_update} \\
        --threads ${threads} \\
        --make-pgen \\
        --out "${study_name}_chr${chr}"
    fi
    """
}
