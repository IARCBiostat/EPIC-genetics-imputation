process PHASING {
    tag "${study_id}_chr${chr}"
    label 'high_mem'
    publishDir "${params.outdir}/${study_id}/phasing", mode: 'copy'

    input:
    tuple val(study_id), val(chr), path(vcf), path(genetic_map)

    output:
    tuple val(study_id), val(chr), path("${study_id}_chr${chr}_phased.vcf.gz"), emit: phased_vcf

    script:
    """
    shapeit4 --input ${vcf} --map ${genetic_map} --region ${chr} --output ${study_id}_chr${chr}_phased.vcf.gz --thread 8
    """
}
