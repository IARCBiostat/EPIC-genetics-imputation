process FILTER_CHROM {
    tag "${study_name}/chr${chr}"

    input:
    tuple val(study_name), val(chr), path(stage2_vcf), path(stage2_tbi)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.filtered.bcf"), path("${study_name}_chr${chr}.id_map.txt")

    script:
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

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
    """
}

process ANNOTATE_CHROM {
    tag "${study_name}/chr${chr}"

    input:
    tuple val(study_name), val(chr), path(filtered_bcf), path(dbsnp_vcf), path(dbsnp_tbi), path(id_map)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.annotated.bcf"), path("${study_name}_chr${chr}.variant_id_map.tsv.gz"), path("${study_name}_chr${chr}.id_map.txt")

    script:
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    set -o pipefail

    \$BCFTOOLS_BIN index "${filtered_bcf}"

    \$BCFTOOLS_BIN annotate \\
      --threads ${threads} \\
      -a "${dbsnp_vcf}" \\
      -c ID \\
      -Ov \\
      -o - \\
      "${filtered_bcf}" \\
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
    def backoff_secs = (task.attempt - 1) * 30
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
        --out ${study_name}_chr${chr}.import
    else
      \$PLINK2_BIN \\
        --bcf ${annotated_bcf} dosage=HDS \\
        --double-id \\
        --update-ids ${id_map} \\
        --threads ${threads} \\
        --make-pgen \\
        --out ${study_name}_chr${chr}.import
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
    tuple val(study_name), val(chr), path(import_pgen), path(import_pvar), path(import_psam), path(id_map)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.pgen"), path("${study_name}_chr${chr}.pvar"), path("${study_name}_chr${chr}.psam"), path("${study_name}_chr${chr}.variant_metrics.tsv"), path("${study_name}_chr${chr}.variant_id_map.tsv.gz")

    script:
    def import_prefix = import_pgen.baseName
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    POST_R2_MAF_VARIANTS=\$(grep -vc '^#' ${import_prefix}.pvar)

    HWE_APPLIED=0
    if [ "${params.run_hwe}" = "true" ] && [ "${chr}" != "X" ]; then
      \$PLINK2_BIN \\
        --pfile ${import_prefix} \\
        --hwe ${params.hwe_p} ${params.hwe_k} midp keep-fewhet \\
        --threads ${threads} \\
        --make-pgen \\
        --out ${study_name}_chr${chr}
      HWE_APPLIED=1
    else
      \$PLINK2_BIN \\
        --pfile ${import_prefix} \\
        --threads ${threads} \\
        --make-pgen \\
        --out ${study_name}_chr${chr}
    fi

    FINAL_VARIANTS=\$(grep -vc '^#' ${study_name}_chr${chr}.pvar)

    read -r INPUT_VARIANTS RSID_COUNT FALLBACK_COUNT DUPLICATE_RSID_FALLBACK_COUNT <<EOF_COUNTS
\$(gzip -cd ${id_map} | awk -F'\\t' '
  NR == 1 { next }
  { input++ }
  \$7 == "rsid" { rsid++ }
  \$7 == "fallback" { fallback++ }
  \$7 == "duplicate_rsid_fallback" { dup++; fallback++ }
  END { printf "%d %d %d %d\\n", input + 0, rsid + 0, fallback + 0, dup + 0 }
')
EOF_COUNTS

    {
      printf "study\\tchr\\tinput_variants\\trsid_variants\\tfallback_variants\\tduplicate_rsid_fallbacks\\tpost_r2_maf_variants\\tfinal_variants\\thwe_applied\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "${chr}" "\${INPUT_VARIANTS}" "\${RSID_COUNT}" "\${FALLBACK_COUNT}" "\${DUPLICATE_RSID_FALLBACK_COUNT}" "\${POST_R2_MAF_VARIANTS}" "\${FINAL_VARIANTS}" "\${HWE_APPLIED}"
    } > ${study_name}_chr${chr}.variant_metrics.tsv
    """
}
