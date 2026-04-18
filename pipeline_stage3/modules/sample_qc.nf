process SAMPLE_QC {
    tag "${study_name}"
    cpus 2
    memory '24 GB'
    publishDir "${params.outdir}/${study_name}/stage3/qc/sample_qc", mode: 'copy', pattern: "${study_name}.*", overwrite: true

    input:
    tuple val(study_name), val(chroms), path(pgens), path(pvars), path(psams), path(variant_stats), path(id_maps), path(stage1_fam)

    output:
    tuple val(study_name), path("${study_name}_allchr_sex.pgen"), path("${study_name}_allchr_sex.pvar"), path("${study_name}_allchr_sex.psam"), path("${study_name}.samples_to_remove.id"), path("${study_name}.sample_qc.tsv"), emit: finalize_input
    path("${study_name}.*"), emit: qc_files

    script:
    def chrom_list = chroms.collect { it.toString() }.join(' ')
    def pgen_list = pgens.collect { it.toString() }.join(' ')
    def pvar_list = pvars.collect { it.toString() }.join(' ')
    def psam_list = psams.collect { it.toString() }.join(' ')
    """
    CHROMS=(${chrom_list})
    PGENS=(${pgen_list})
    PVARS=(${pvar_list})
    PSAMS=(${psam_list})

    : > merge_list.txt
    for idx in "\${!PGENS[@]}"; do
      printf "%s %s %s\\n" "\${PGENS[\$idx]}" "\${PVARS[\$idx]}" "\${PSAMS[\$idx]}" >> merge_list.txt
    done

    \$PLINK2_BIN \\
      --pmerge-list merge_list.txt \\
      --make-pgen \\
      --sort-vars \\
      --out ${study_name}_allchr

    awk '{print \$1, \$2, \$5}' "${stage1_fam}" > sex_update.txt

    \$PLINK2_BIN \\
      --pfile ${study_name}_allchr \\
      --update-sex sex_update.txt \\
      --make-pgen \\
      --out ${study_name}_allchr_sex

    \$PLINK2_BIN --pfile ${study_name}_allchr_sex --make-bed --out ${study_name}_allchr_bed
    \$PLINK2_BIN --pfile ${study_name}_allchr_sex --chr 1-22 --make-bed --out ${study_name}_autosomes_bed

    \$PLINK_BIN --bfile ${study_name}_autosomes_bed --indep-pairwise 1500 150 0.2 --out ${study_name}_pruned
    \$PLINK_BIN --bfile ${study_name}_autosomes_bed --extract ${study_name}_pruned.prune.in --make-bed --out ${study_name}_pruned_bed

    if ! \$PLINK_BIN --bfile ${study_name}_allchr_bed --check-sex --out ${study_name}_sexcheck; then
      : > ${study_name}_sexcheck.sexcheck
    fi

    if ! \$PLINK2_BIN --bfile ${study_name}_pruned_bed --make-king-table --out ${study_name}_king; then
      : > ${study_name}_related.king.cutoff.out.id
    fi

    if [ -f ${study_name}_king.kin0 ]; then
      if ! \$PLINK2_BIN --king-cutoff-table ${study_name}_king.kin0 ${params.king_cutoff} --out ${study_name}_related; then
        : > ${study_name}_related.king.cutoff.out.id
      fi
    else
      : > ${study_name}_related.king.cutoff.out.id
    fi

    \$PLINK_BIN --bfile ${study_name}_pruned_bed --het --out ${study_name}_het
    \$PLINK2_BIN --bfile ${study_name}_pruned_bed --pca ${params.ancestry_pc_count} --out ${study_name}_pca

    \$PYTHON3_BIN "${projectDir}/bin/identify_sample_outliers.py" \\
      --eigenvec ${study_name}_pca.eigenvec \\
      --het ${study_name}_het.het \\
      --sexcheck ${study_name}_sexcheck.sexcheck \\
      --out-prefix ${study_name} \\
      --pc-count ${params.ancestry_pc_count} \\
      --ancestry-z-threshold ${params.ancestry_z_threshold} \\
      --het-sd-threshold ${params.het_sd_threshold}

    if [ -f ${study_name}_related.king.cutoff.out.id ]; then
      cp ${study_name}_related.king.cutoff.out.id ${study_name}.related_outliers.id
    else
      : > ${study_name}.related_outliers.id
    fi

    : > ${study_name}.samples_to_remove.id
    cat ${study_name}.sex_mismatch.id >> ${study_name}.samples_to_remove.id
    cat ${study_name}.heterozygosity_outliers.id >> ${study_name}.samples_to_remove.id
    cat ${study_name}.related_outliers.id >> ${study_name}.samples_to_remove.id

    ANCESTRY_REMOVED_COUNT=0
    if [ "${params.exclude_ancestry_outliers}" = "true" ]; then
      cat ${study_name}.ancestry_outliers.id >> ${study_name}.samples_to_remove.id
      ANCESTRY_REMOVED_COUNT=\$(wc -l < ${study_name}.ancestry_outliers.id | tr -d ' ')
    fi

    sort -u ${study_name}.samples_to_remove.id -o ${study_name}.samples_to_remove.id

    PRE_FINAL_SAMPLES=\$(awk 'NR > 1 {count++} END {print count + 0}' ${study_name}_allchr_sex.psam)
    SEX_MISMATCH_COUNT=\$(wc -l < ${study_name}.sex_mismatch.id | tr -d ' ')
    RELATED_COUNT=\$(wc -l < ${study_name}.related_outliers.id | tr -d ' ')
    HET_COUNT=\$(wc -l < ${study_name}.heterozygosity_outliers.id | tr -d ' ')
    ANCESTRY_IDENTIFIED_COUNT=\$(wc -l < ${study_name}.ancestry_outliers.id | tr -d ' ')
    TOTAL_REMOVED_COUNT=\$(wc -l < ${study_name}.samples_to_remove.id | tr -d ' ')

    {
      printf "study\\tpre_final_samples\\tsex_mismatch\\trelated_removed\\theterozygosity_outliers\\tancestry_outliers_identified\\tancestry_outliers_removed\\ttotal_removed\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "\${PRE_FINAL_SAMPLES}" "\${SEX_MISMATCH_COUNT}" "\${RELATED_COUNT}" "\${HET_COUNT}" "\${ANCESTRY_IDENTIFIED_COUNT}" "\${ANCESTRY_REMOVED_COUNT}" "\${TOTAL_REMOVED_COUNT}"
    } > ${study_name}.sample_qc.tsv
    """
}
