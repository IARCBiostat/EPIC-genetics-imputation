// pipeline/modules/imputation.nf

process IMPUTE_AUTOSOMES {
    tag "${study_name}/chr${chr}"
    cpus { params.minimac_threads }
    memory '32 GB'
    publishDir "${params.outdir}/${study_name}/stage2", mode: 'copy', overwrite: true

    input:
    tuple val(study_name), val(chr), path(phased_vcf), path(phased_vcf_index), path(msav_panel)
    
    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}_GxS.imputed.vcf.gz"), path("${study_name}_chr${chr}_GxS.imputed.vcf.gz.tbi")

    script:
    """
    \$MINIMAC4_BIN \\
        ${msav_panel} \\
        ${phased_vcf} \\
        -o ${study_name}_chr${chr}_GxS.imputed.vcf.gz \\
        --threads ${params.minimac_threads} \\
        -b ${params.minimac_batch_size} \\
        --min-r2 ${params.min_r2} \\
        -O vcf.gz
        
    \$BCFTOOLS_BIN index -t ${study_name}_chr${chr}_GxS.imputed.vcf.gz
    """
}
