process PREPARE_CHROM {
    tag "${study_name}/chr${chr}"

    input:
    tuple val(study_name), val(chr), path(stage2_vcf), path(stage2_tbi), path(dbsnp_vcf), path(dbsnp_tbi)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.annotated.bcf"), path("${study_name}_chr${chr}.variant_id_map.tsv.gz"), path("${study_name}_chr${chr}.id_map.txt")

    script:
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    set -o pipefail

    \$BCFTOOLS_BIN query -l "${stage2_vcf}" | awk '{
        full_id = \$1;
        n = length(full_id);
        half = (n - 1) / 2;
        left = substr(full_id, 1, half);
        right = substr(full_id, half + 2);
        if (left == right) {
            print full_id, full_id, left, left;
        } else {
            print full_id, full_id, full_id, full_id;
        }
    }' > ${study_name}_chr${chr}.id_map.txt

    \$BCFTOOLS_BIN view \\
      --threads ${threads} \\
      -i 'INFO/R2>=${params.min_r2} && INFO/MAF>=${params.maf}' \\
      -Ob \\
      -o ${study_name}_chr${chr}.filtered.bcf \\
      "${stage2_vcf}"

    \$BCFTOOLS_BIN index "${study_name}_chr${chr}.filtered.bcf"

    \$BCFTOOLS_BIN annotate \\
      --threads ${threads} \\
      -a "${dbsnp_vcf}" \\
      -c ID \\
      -Ov \\
      -o - \\
      "${study_name}_chr${chr}.filtered.bcf" \\
    | \$PYTHON3_BIN "${projectDir}/bin/normalize_variant_ids.py" \\
      --input - \\
      --output - \\
      --mapping-output ${study_name}_chr${chr}.variant_id_map.tsv.gz \\
    | \$BCFTOOLS_BIN view \\
      --threads ${threads} \\
      -Ob \\
      -o ${study_name}_chr${chr}.annotated.bcf \\
      -
    """
}

process IMPORT_CHROM {
    tag "${study_name}/chr${chr}"

    input:
    tuple val(study_name), val(chr), path(annotated_bcf), path(variant_id_map), path(id_map)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.import.pgen"), path("${study_name}_chr${chr}.import.pvar"), path("${study_name}_chr${chr}.import.psam"), path("${study_name}_chr${chr}.variant_id_map.tsv.gz")

    script:
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 120
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    if [ "${chr}" = "X" ]; then
      \$PLINK2_BIN \\
        --bcf ${annotated_bcf} dosage=HDS \\
        --double-id \\
        --split-par b38 \\
        --lax-chrx-import \\
        --update-ids ${id_map} \\
        --threads ${threads} \\
        --make-pgen \\
        --out "${study_name}_chr${chr}.import"
    else
      \$PLINK2_BIN \\
        --bcf ${annotated_bcf} dosage=HDS \\
        --double-id \\
        --update-ids ${id_map} \\
        --threads ${threads} \\
        --make-pgen \\
        --out "${study_name}_chr${chr}.import"
    fi
    """
}

process HWE_CHROM {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith('.variant_metrics.tsv')) {
            return "${study_name}/stage3/report/tables/${filename}"
        }
        if (filename.endsWith('.variant_id_map.tsv.gz')) {
            return "${study_name}/stage3/report/manifests/${filename}"
        }
        if (params.publish_intermediate_plink.toString().toBoolean() && (filename.endsWith('.pgen') || filename.endsWith('.pvar') || filename.endsWith('.psam'))) {
            return "${study_name}/stage3/prep_chrom/${filename}"
        }
        null
    }

    input:
    tuple val(study_name), val(chr), path(import_pgen), path(import_pvar), path(import_psam), path(variant_id_map), path(stage1_fam)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.pgen"), path("${study_name}_chr${chr}.pvar"), path("${study_name}_chr${chr}.psam"), path("${study_name}_chr${chr}.variant_metrics.tsv"), emit: chrom_data
    path("${study_name}_chr${chr}.variant_id_map.tsv.gz")
    tuple val(study_name), path("${study_name}_chr${chr}.hwe.exclude"), emit: hwe_exclude

    script:
    def import_prefix = import_pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    POST_R2_MAF_VARIANTS=\$(grep -vc '^#' ${import_prefix}.pvar)

    HWE_APPLIED=0
    HWE_NO_CONTROLS=0
    HWE_EXCLUDE_COUNT=0

    if [ "${params.hwe}" = "true" ] && [ "${chr}" != "X" ]; then
      # Build pheno.txt using the PSAM's FID/IID (= the PGEN's sample IDs) so that
      # PLINK2 can match samples.  PLINK1 --recode vcf emits "FID_IID" for samples
      # where FID != IID, so after --double-id the PGEN has FID=IID="FID_IID".
      # We index the FAM by IID alone and by the "FID_IID" concatenation to cover
      # both cases.
      awk '
        NR==FNR {
          pheno_iid[\$2] = \$6
          pheno_fid_iid[\$1 "_" \$2] = \$6
          next
        }
        /^#/ { next }
        {
          iid = \$2
          if (iid in pheno_iid)        p = pheno_iid[iid]
          else if (iid in pheno_fid_iid) p = pheno_fid_iid[iid]
          else                           p = -9
          print \$1, \$2, p
        }
      ' "${stage1_fam}" "${import_prefix}.psam" > pheno.txt
      awk '\$3 == 1 {print \$1, \$2}' pheno.txt > controls.txt
      CONTROL_COUNT=\$(wc -l < controls.txt | tr -d ' ')
      if [ "\${CONTROL_COUNT}" -ge 2 ]; then
        HWE_STEP_FAILED=0
        \$PLINK2_BIN \\
          --pfile ${import_prefix} \\
          --keep controls.txt \\
          --hwe ${params.hwe_p} midp keep-fewhet \\
          --write-snplist allow-dups \\
          --threads ${threads} \\
          --out hwe_pass || HWE_STEP_FAILED=1

        if [ "\${HWE_STEP_FAILED}" -eq 0 ] && [ -s hwe_pass.snplist ]; then
          awk '!/^#/{print \$3}' ${import_prefix}.pvar | sort > _all_variants.txt
          sort hwe_pass.snplist > _hwe_pass_sorted.txt
          comm -23 _all_variants.txt _hwe_pass_sorted.txt > ${study_name}_chr${chr}.hwe.exclude
        else
          echo "WARNING: HWE step failed or produced empty snplist for ${study_name} chr${chr}; no variants excluded." >&2
          : > ${study_name}_chr${chr}.hwe.exclude
        fi
        HWE_EXCLUDE_COUNT=\$(wc -l < ${study_name}_chr${chr}.hwe.exclude | tr -d ' ')
        HWE_APPLIED=1
      else
        echo "WARNING: HWE filtering requested but fewer than 2 control samples found for ${study_name} chr${chr} (n=\${CONTROL_COUNT}); no exclusion list generated." >&2
        HWE_NO_CONTROLS=1
        : > ${study_name}_chr${chr}.hwe.exclude
      fi
    else
      : > ${study_name}_chr${chr}.hwe.exclude
    fi

    \$PLINK2_BIN \\
      --pfile ${import_prefix} \\
      --threads ${threads} \\
      --make-pgen \\
      --out ${study_name}_chr${chr}

    FINAL_VARIANTS=\$(grep -vc '^#' ${study_name}_chr${chr}.pvar)

    read -r INPUT_VARIANTS RSID_COUNT FALLBACK_COUNT DUPLICATE_RSID_FALLBACK_COUNT <<EOF_COUNTS
\$(gzip -cd ${variant_id_map} | awk -F'\\t' '
  NR == 1 { next }
  { input++ }
  \$7 == "rsid" { rsid++ }
  \$7 == "fallback" { fallback++ }
  \$7 == "duplicate_rsid_fallback" { dup++; fallback++ }
  END { printf "%d %d %d %d\\n", input + 0, rsid + 0, fallback + 0, dup + 0 }
')
EOF_COUNTS

    {
      printf "study\\tchr\\tinput_variants\\trsid_variants\\tfallback_variants\\tduplicate_rsid_fallbacks\\tpost_r2_maf_variants\\tfinal_variants\\thwe_applied\\thwe_no_controls\\thwe_exclude_count\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "${chr}" "\${INPUT_VARIANTS}" "\${RSID_COUNT}" "\${FALLBACK_COUNT}" "\${DUPLICATE_RSID_FALLBACK_COUNT}" "\${POST_R2_MAF_VARIANTS}" "\${FINAL_VARIANTS}" "\${HWE_APPLIED}" "\${HWE_NO_CONTROLS}" "\${HWE_EXCLUDE_COUNT}"
    } > ${study_name}_chr${chr}.variant_metrics.tsv

    test ${variant_id_map} -ef ${study_name}_chr${chr}.variant_id_map.tsv.gz || cp ${variant_id_map} ${study_name}_chr${chr}.variant_id_map.tsv.gz
    """
}

process AGGREGATE_HWE_EXCLUDE {
    tag "${study_name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${study_name}/stage3/report/flags/${filename}"
    }

    input:
    tuple val(study_name), path(hwe_excludes)

    output:
    tuple val(study_name), path("${study_name}.hwe.exclude")

    script:
    """
    cat ${hwe_excludes} | sort -u > ${study_name}.hwe.exclude
    """
}
