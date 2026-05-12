// pipeline/modules/phasing_metrics.nf

process STAGE_PHASING_METRICS {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        return "${study_name}/stage2/report/tables/phasing/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(phased_vcf), path(phased_vcf_index)

    output:
    tuple path("chr${chr}.phase_quality.tsv"), path("chr${chr}.phasing_metrics.tsv")

    script:
    """
    \$PYTHON3_BIN ${projectDir}/bin/write_phasing_metrics.py \\
        --vcf ${phased_vcf} \\
        --study ${study_name} \\
        --chrom ${chr} \\
        --variant-output chr${chr}.phase_quality.tsv \\
        --summary-output chr${chr}.phasing_metrics.tsv
    """
}
