process PRUNE_CHROM_FOR_QC {
    label 'sample_review'
    tag "${study_name}/chr${chr}"

    input:
    tuple val(study_name), val(chr), path(pgen), path(pvar), path(psam)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}_pruned.bed"), path("${study_name}_chr${chr}_pruned.bim"), path("${study_name}_chr${chr}_pruned.fam"), emit: pruned_chr

    script:
    def pfile_prefix = pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    \$PLINK2_BIN \\
      --pfile ${pfile_prefix} \\
      --set-all-var-ids '@:#:\$r:\$a' \\
      --new-id-max-allele-len 1000 missing \\
      --rm-dup force-first \\
      --indep-pairwise 1500 150 0.2 \\
      --threads ${threads} \\
      --out ${study_name}_chr${chr}_prune

    \$PLINK2_BIN \\
      --pfile ${pfile_prefix} \\
      --set-all-var-ids '@:#:\$r:\$a' \\
      --new-id-max-allele-len 1000 missing \\
      --rm-dup force-first \\
      --extract ${study_name}_chr${chr}_prune.prune.in \\
      --threads ${threads} \\
      --make-bed \\
      --out ${study_name}_chr${chr}_pruned
    """
}

process MERGE_PRUNED_FOR_QC {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (params.publish_intermediate_plink.toString().toBoolean() && (filename.endsWith('.bed') || filename.endsWith('.bim') || filename.endsWith('.fam'))) {
            return "${study_name}/stage3/sample_review/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(beds), path(bims), path(fams)

    output:
    tuple val(study_name), path("${study_name}_pruned_bed.bed"), path("${study_name}_pruned_bed.bim"), path("${study_name}_pruned_bed.fam"), emit: pruned_bed

    script:
    def bed_prefixes = beds.collect { it.baseName }.join(' ')
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    : > merge_list.txt
    for pfx in ${bed_prefixes}; do
      echo "\${pfx}" >> merge_list.txt
    done

    \$PLINK2_BIN \\
      --pmerge-list merge_list.txt bfile \\
      --threads ${threads} \\
      --memory 20000 \\
      --make-bed \\
      --out "${study_name}_pruned_bed"
    """
}

process MAKE_UPDATE_FILES {
    label 'sample_review'
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.sex_update.txt') || filename.endsWith('.pheno_update.txt')) {
            return "${study_name}/stage3/report/manifests/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(chr1_psam), path(stage1_fam)

    output:
    tuple val(study_name), path("${study_name}.sex_update.txt"), path("${study_name}.pheno_update.txt"), emit: update_files

    script:
    """
    awk '
      NR==FNR {
        sex_by_iid[\$2] = \$5
        sex_by_fid_iid[\$1 "_" \$2] = \$5
        next
      }
      /^#/ { next }
      {
        iid = \$2
        if (iid in sex_by_iid)          s = sex_by_iid[iid]
        else if (iid in sex_by_fid_iid) s = sex_by_fid_iid[iid]
        else                             s = 0
        print \$1, \$2, s
      }
    ' "${stage1_fam}" ${chr1_psam} > ${study_name}.sex_update.txt

    awk '
      BEGIN { print "#FID\tIID\tPHENO1" }
      NR==FNR {
        pheno_by_iid[\$2] = \$6
        pheno_by_fid_iid[\$1 "_" \$2] = \$6
        next
      }
      /^#/ { next }
      {
        iid = \$2
        if (iid in pheno_by_iid)          p = pheno_by_iid[iid]
        else if (iid in pheno_by_fid_iid) p = pheno_by_fid_iid[iid]
        else                               p = -9
        print \$1, \$2, p
      }
    ' "${stage1_fam}" ${chr1_psam} >> ${study_name}.pheno_update.txt
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
        if (filename.endsWith('.id') || filename.endsWith('.exclude')) {
            return "${study_name}/stage3/report/flags/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), path(related_ids), path(het), path(eigenvec), path(eigenval), path(chr1_psam)
    val exclude_ancestry_outliers
    val exclude_related

    output:
    tuple val(study_name), path("${study_name}.samples_to_remove.id"), path("${study_name}.sample_review.tsv"), emit: summary
    tuple val(study_name), path("${study_name}.related.exclude"), emit: related_exclude
    tuple val(study_name), path("${study_name}.ancestry.exclude"), emit: ancestry_exclude
    path("${study_name}.heterozygosity_outliers.id"), emit: heterozygosity_outliers

    script:
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    \$PYTHON3_BIN "${projectDir}/bin/identify_sample_outliers.py" \\
      --eigenvec ${eigenvec} \\
      --het ${het} \\
      --out-prefix ${study_name} \\
      --pc-count ${params.ancestry_pc_count} \\
      --ancestry-z-threshold ${params.ancestry_z_threshold} \\
      --het-sd-threshold ${params.het_sd_threshold}

    if [ -f ${related_ids} ]; then
      cp ${related_ids} ${study_name}.related.exclude
    else
      : > ${study_name}.related.exclude
    fi

    : > ${study_name}.samples_to_remove.id
    cat ${study_name}.heterozygosity_outliers.id >> ${study_name}.samples_to_remove.id

    RELATED_REMOVED_COUNT=0
    if [ "${exclude_related}" = "true" ]; then
      cat ${study_name}.related.exclude >> ${study_name}.samples_to_remove.id
      RELATED_REMOVED_COUNT=\$(wc -l < ${study_name}.related.exclude | tr -d ' ')
    fi

    cp ${study_name}.ancestry_outliers.id ${study_name}.ancestry.exclude

    ANCESTRY_REMOVED_COUNT=0
    if [ "${exclude_ancestry_outliers}" = "true" ]; then
      cat ${study_name}.ancestry.exclude >> ${study_name}.samples_to_remove.id
      ANCESTRY_REMOVED_COUNT=\$(wc -l < ${study_name}.ancestry.exclude | tr -d ' ')
    fi

    sort -u ${study_name}.samples_to_remove.id -o ${study_name}.samples_to_remove.id

    PRE_FINAL_SAMPLES=\$(awk 'NR > 1 {count++} END {print count + 0}' ${chr1_psam})
    RELATED_IDENTIFIED_COUNT=\$(wc -l < ${study_name}.related.exclude | tr -d ' ')
    HET_COUNT=\$(wc -l < ${study_name}.heterozygosity_outliers.id | tr -d ' ')
    ANCESTRY_IDENTIFIED_COUNT=\$(wc -l < ${study_name}.ancestry.exclude | tr -d ' ')
    TOTAL_REMOVED_COUNT=\$(wc -l < ${study_name}.samples_to_remove.id | tr -d ' ')

    {
      printf "study\\tpre_final_samples\\trelated_identified\\trelated_removed\\theterozygosity_outliers\\tancestry_outliers_identified\\tancestry_outliers_removed\\ttotal_removed\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "\${PRE_FINAL_SAMPLES}" "\${RELATED_IDENTIFIED_COUNT}" "\${RELATED_REMOVED_COUNT}" "\${HET_COUNT}" "\${ANCESTRY_IDENTIFIED_COUNT}" "\${ANCESTRY_REMOVED_COUNT}" "\${TOTAL_REMOVED_COUNT}"
    } > ${study_name}.sample_review.tsv
    """
}
