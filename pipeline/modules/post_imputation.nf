process POST_IMPUTATION_VARIANT_QC {
    tag "${study_id}_chr${chr}"
    label 'mid_mem'
    publishDir "${params.outdir}/${study_id}/post_imputation/variants", mode: 'copy'

    input:
    tuple val(study_id), val(chr), path(imputed_vcf)

    output:
    tuple val(study_id), val(chr), path("${study_id}_chr${chr}_filtered.vcf.gz"), emit: filtered_vcf
    tuple val(study_id), val(chr), path("${study_id}_chr${chr}_filtered.bed"), path("${study_id}_chr${chr}_filtered.bim"), path("${study_id}_chr${chr}_filtered.fam"), emit: filtered_plink

    script:
    """
    # 1. Filter by R2 > 0.7 and MAF > 0.01
    bcftools filter -i 'R2>${params.info_score} && MAF>${params.maf}' ${imputed_vcf} -Oz -o ${study_id}_chr${chr}_bcftools_filtered.vcf.gz
    
    # 2. Convert to PLINK and apply HWE
    plink --vcf ${study_id}_chr${chr}_bcftools_filtered.vcf.gz \
          --hwe ${params.hwe} \
          --make-bed \
          --out ${study_id}_chr${chr}_filtered
          
    # Convert back to VCF for consistency
    plink --bfile ${study_id}_chr${chr}_filtered --recode vcf --out ${study_id}_chr${chr}_filtered
    bgzip ${study_id}_chr${chr}_filtered.vcf
    """
}

process POST_IMPUTATION_SAMPLE_QC {
    tag "${study_id}"
    label 'high_mem'
    publishDir "${params.outdir}/${study_id}/post_imputation/samples", mode: 'copy'

    input:
    tuple val(study_id), path(beds), path(bims), path(fams)

    output:
    path "${study_id}_qc_pass_samples.txt", emit: pass_samples
    path "${study_id}_sample_qc_report.txt", emit: report

    script:
    """
    # 1. Merge all chromosomes
    echo "${bims}" | tr ' ' '\n' | sed 's/\.bim//' > merge_list.txt
    plink --merge-list merge_list.txt --make-bed --out ${study_id}_merged
    
    # 2. LD Pruning for sample QC
    plink --bfile ${study_id}_merged --indep-pairwise 1500 150 0.2 --out pruned
    plink --bfile ${study_id}_merged --extract pruned.prune.in --make-bed --out ${study_id}_pruned
    
    # 3. Sex Check
    plink --bfile ${study_id}_pruned --check-sex --out ${study_id}_sex
    
    # 4. Relatedness (KING)
    plink --bfile ${study_id}_pruned --make-king-table --out ${study_id}_king
    plink --bfile ${study_id}_pruned --king-cutoff 0.088 --out ${study_id}_related
    
    # 5. Heterozygosity
    plink --bfile ${study_id}_pruned --het --out ${study_id}_het
    
    # Placeholder for outlier removal logic
    cat ${study_id}_merged.fam | awk '{print $1, $2}' > ${study_id}_qc_pass_samples.txt
    echo "Post-imputation sample QC complete." > ${study_id}_sample_qc_report.txt
    """
}
