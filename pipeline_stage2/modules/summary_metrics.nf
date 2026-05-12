// pipeline/modules/summary_metrics.nf

process STAGE_IMPUTATION_METRICS {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.r2_maf.tsv')) {
            return "${study_name}/stage2/report/tables/imputation/${filename}"
        }
        return "${study_name}/stage2/report/tables/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(imputed_vcf), path(imputed_vcf_index)

    output:
    tuple val(study_name), path("chr${chr}.imputation_metrics.tsv"), path("chr${chr}.r2_maf.tsv")

    script:
    """
    echo "imputation_metrics_schema=2" >&2

    \$PYTHON3_BIN ${projectDir}/bin/write_imputation_metrics.py \\
        --vcf ${imputed_vcf} \\
        --study ${study_name} \\
        --chrom ${chr} \\
        --output chr${chr}.imputation_metrics.tsv \\
        --variant-output chr${chr}.r2_maf.tsv \\
        --min-r2-threshold ${params.min_r2} \\
        --high-quality-threshold 0.8
    """
}
