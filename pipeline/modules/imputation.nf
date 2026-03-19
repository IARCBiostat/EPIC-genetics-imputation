process IMPUTATION {
    tag "${study_id}_chr${chr}"
    label 'high_mem'
    publishDir "${params.outdir}/${study_id}/imputation", mode: 'copy'

    input:
    tuple val(study_id), val(chr), path(phased_vcf), path(ref_vcf)

    output:
    tuple val(study_id), val(chr), path("${study_id}_chr${chr}_imputed.vcf.gz"), emit: imputed_vcf

    script:
    """
    minimac4 --refHaps ${ref_vcf} --haps ${phased_vcf} --prefix ${study_id}_chr${chr}_imputed --format GT,DS,GP
    """
}
