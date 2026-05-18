process MERGE_STUDY {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.sex_update.txt')) {
            return "${study_name}/stage3/report/manifests/${filename}"
        }
        if (params.publish_intermediate_plink.toString().toBoolean() && (filename.endsWith('.pgen') || filename.endsWith('.pvar') || filename.endsWith('.psam'))) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), val(chroms), path(pgens), path(pvars), path(psams), path(variant_stats), path(id_maps), path(stage1_fam)

    output:
    tuple val(study_name), val(chroms), path("${study_name}_allchr.pgen"), path("${study_name}_allchr.pvar"), path("${study_name}_allchr.psam"), path("${study_name}.sex_update.txt"), emit: merged

    script:
    def pgen_list = pgens.collect { it.toString() }.join(' ')
    def pvar_list = pvars.collect { it.toString() }.join(' ')
    def psam_list = psams.collect { it.toString() }.join(' ')
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    PGENS=(${pgen_list})
    PVARS=(${pvar_list})
    PSAMS=(${psam_list})

    : > merge_list.txt
    for idx in "\${!PGENS[@]}"; do
      printf "%s %s %s\\n" "\${PGENS[\$idx]}" "\${PVARS[\$idx]}" "\${PSAMS[\$idx]}" >> merge_list.txt
    done

    \$PLINK2_BIN \\
      --pmerge-list merge_list.txt \\
      --threads ${threads} \\
      --make-pgen \\
      --sort-vars \\
      --out ${study_name}_allchr

    awk '{print \$1, \$2, \$5}' "${stage1_fam}" > ${study_name}.sex_update.txt
    """
}

process SEX_CHECK {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.sexcheck')) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        if (params.publish_intermediate_plink.toString().toBoolean() && (filename.endsWith('.bed') || filename.endsWith('.bim') || filename.endsWith('.fam'))) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), val(chroms), path(merged_pgen), path(merged_pvar), path(merged_psam), path(sex_update)

    output:
    tuple val(study_name), path("${study_name}_sexcheck.sexcheck"), emit: sexcheck

    script:
    def chrom_list = chroms.collect { it.toString() }.join(' ')
    def merged_prefix = merged_pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    CHROMS=(${chrom_list})

    HAS_CHRX=0
    for chr in "\${CHROMS[@]}"; do
      if [ "\${chr}" = "X" ]; then
        HAS_CHRX=1
      fi
    done

    if [ "\${HAS_CHRX}" = "1" ]; then
      if ! \$PLINK2_BIN \\
          --pfile ${merged_prefix} \\
          --chr X \\
          --update-sex ${sex_update} \\
          --threads ${threads} \\
          --make-bed \\
          --out ${study_name}_chrX_bed; then
        : > ${study_name}_sexcheck.sexcheck
      elif ! \$PLINK_BIN --bfile ${study_name}_chrX_bed --check-sex --out ${study_name}_sexcheck; then
        : > ${study_name}_sexcheck.sexcheck
      fi
    else
      : > ${study_name}_sexcheck.sexcheck
    fi
    """
}

process PRUNE_AUTOSOMES {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (params.publish_intermediate_plink.toString().toBoolean() && (filename.endsWith('.bed') || filename.endsWith('.bim') || filename.endsWith('.fam'))) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(merged_pgen), path(merged_pvar), path(merged_psam)

    output:
    tuple val(study_name), path("${study_name}_pruned_bed.bed"), path("${study_name}_pruned_bed.bim"), path("${study_name}_pruned_bed.fam"), emit: pruned_bed

    script:
    def merged_prefix = merged_pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    \$PLINK2_BIN \\
      --pfile ${merged_prefix} \\
      --chr 1-22 \\
      --indep-pairwise 1500 150 0.2 \\
      --threads ${threads} \\
      --out ${study_name}_pruned

    \$PLINK2_BIN \\
      --pfile ${merged_prefix} \\
      --chr 1-22 \\
      --extract ${study_name}_pruned.prune.in \\
      --threads ${threads} \\
      --make-bed \\
      --out ${study_name}_pruned_bed
    """
}

process KING_QC {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.kin0')) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        if (filename.endsWith('.id')) {
            return "${study_name}/stage3/report/flags/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(pruned_bed), path(pruned_bim), path(pruned_fam)

    output:
    tuple val(study_name), path("${study_name}_related.king.cutoff.out.id"), emit: related_ids
    path("${study_name}_king.kin0"), emit: king_table

    script:
    def pruned_prefix = pruned_bed.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    if ! \$PLINK2_BIN --bfile ${pruned_prefix} --threads ${threads} --make-king-table --out ${study_name}_king; then
      : > ${study_name}_king.kin0
      : > ${study_name}_related.king.cutoff.out.id
    fi

    [ -f ${study_name}_king.kin0 ] || : > ${study_name}_king.kin0

    if [ -s ${study_name}_king.kin0 ]; then
      if ! \$PLINK2_BIN --king-cutoff-table ${study_name}_king.kin0 ${params.king_cutoff} --out ${study_name}_related; then
        : > ${study_name}_related.king.cutoff.out.id
      fi
    else
      : > ${study_name}_related.king.cutoff.out.id
    fi

    [ -f ${study_name}_related.king.cutoff.out.id ] || : > ${study_name}_related.king.cutoff.out.id
    """
}

process HET_PCA_QC {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.het') || filename.endsWith('.eigenvec') || filename.endsWith('.eigenval')) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(pruned_bed), path(pruned_bim), path(pruned_fam)

    output:
    tuple val(study_name), path("${study_name}_het.het"), path("${study_name}_pca.eigenvec"), path("${study_name}_pca.eigenval"), emit: qc_files

    script:
    def pruned_prefix = pruned_bed.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    \$PLINK_BIN --bfile ${pruned_prefix} --het --out ${study_name}_het
    \$PLINK2_BIN --bfile ${pruned_prefix} --threads ${threads} --pca ${params.ancestry_pc_count} --out ${study_name}_pca
    """
}

process SAMPLE_REVIEW_SUMMARY {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.sample_review.tsv')) {
            return "${study_name}/stage3/report/tables/${filename}"
        }
        if (filename.endsWith('.samples_to_remove.id')) {
            return "${study_name}/stage3/report/manifests/${filename}"
        }
        if (filename.endsWith('.id')) {
            return "${study_name}/stage3/report/flags/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(sexcheck), path(related_ids), path(het), path(eigenvec), path(eigenval), path(merged_psam)

    output:
    tuple val(study_name), path("${study_name}.samples_to_remove.id"), path("${study_name}.sample_review.tsv"), emit: summary
    path("${study_name}.sex_mismatch.id"), emit: sex_mismatch
    path("${study_name}.heterozygosity_outliers.id"), emit: heterozygosity_outliers
    path("${study_name}.ancestry_outliers.id"), emit: ancestry_outliers
    path("${study_name}.related_outliers.id"), emit: related_outliers

    script:
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    \$PYTHON3_BIN "${projectDir}/bin/identify_sample_outliers.py" \\
      --eigenvec ${eigenvec} \\
      --het ${het} \\
      --sexcheck ${sexcheck} \\
      --out-prefix ${study_name} \\
      --pc-count ${params.ancestry_pc_count} \\
      --ancestry-z-threshold ${params.ancestry_z_threshold} \\
      --het-sd-threshold ${params.het_sd_threshold}

    if [ -f ${related_ids} ]; then
      cp ${related_ids} ${study_name}.related_outliers.id
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

    PRE_FINAL_SAMPLES=\$(awk 'NR > 1 {count++} END {print count + 0}' ${merged_psam})
    SEX_MISMATCH_COUNT=\$(wc -l < ${study_name}.sex_mismatch.id | tr -d ' ')
    RELATED_COUNT=\$(wc -l < ${study_name}.related_outliers.id | tr -d ' ')
    HET_COUNT=\$(wc -l < ${study_name}.heterozygosity_outliers.id | tr -d ' ')
    ANCESTRY_IDENTIFIED_COUNT=\$(wc -l < ${study_name}.ancestry_outliers.id | tr -d ' ')
    TOTAL_REMOVED_COUNT=\$(wc -l < ${study_name}.samples_to_remove.id | tr -d ' ')

    {
      printf "study\\tpre_final_samples\\tsex_mismatch\\trelated_removed\\theterozygosity_outliers\\tancestry_outliers_identified\\tancestry_outliers_removed\\ttotal_removed\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "\${PRE_FINAL_SAMPLES}" "\${SEX_MISMATCH_COUNT}" "\${RELATED_COUNT}" "\${HET_COUNT}" "\${ANCESTRY_IDENTIFIED_COUNT}" "\${ANCESTRY_REMOVED_COUNT}" "\${TOTAL_REMOVED_COUNT}"
    } > ${study_name}.sample_review.tsv
    """
}
