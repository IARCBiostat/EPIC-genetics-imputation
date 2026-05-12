process STAGE_EMPIRICAL_VALIDATION_METRICS {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        return "${study_name}/stage2/report/tables/imputation_validation/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(target_vcf), path(target_vcf_index), path(imputed_vcf), path(imputed_vcf_index)

    output:
    tuple val(study_name), path("chr${chr}.empirical_metrics.tsv"), path("chr${chr}.empirical_summary.tsv")

    script:
    """
    \$PYTHON3_BIN ${projectDir}/bin/write_empirical_validation_metrics.py \\
        --target-vcf ${target_vcf} \\
        --imputed-vcf ${imputed_vcf} \\
        --study ${study_name} \\
        --chrom ${chr} \\
        --variant-output chr${chr}.empirical_metrics.tsv \\
        --summary-output chr${chr}.empirical_summary.tsv \\
        --min-samples ${params.empirical_min_samples} \\
        --min-dose0-samples ${params.dose0_min_samples}
    """

    stub:
    """
    printf "chrom\\tpos\\tvariant_id\\tmaf\\tn_samples\\tempirical_dosage_r2\\tdose0\\tdose0_n\\tmean_observed_dosage\\tmean_imputed_dosage\\tdosage_bias\\timputation_r2\\n" > chr${chr}.empirical_metrics.tsv
    printf "study\\tchrom\\tvalidation_variants\\tempirical_r2_variants\\tmean_empirical_dosage_r2\\tmedian_empirical_dosage_r2\\tp5_empirical_dosage_r2\\tdose0_variants\\tmean_dose0\\tmedian_dose0\\tmean_abs_dosage_bias\\n${study_name}\\t${chr}\\t0\\t0\\t\\t\\t\\t0\\t\\t\\t\\n" > chr${chr}.empirical_summary.tsv
    """
}
