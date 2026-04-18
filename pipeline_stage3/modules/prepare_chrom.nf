process PREP_CHROM {
    tag "${study_name}/chr${chr}"
    cpus 2
    memory '16 GB'
    publishDir "${params.outdir}/${study_name}/stage3/qc/variant_qc", mode: 'copy', pattern: "${study_name}_chr${chr}.variant_qc.tsv", overwrite: true
    publishDir "${params.outdir}/${study_name}/stage3/qc/variant_qc", mode: 'copy', pattern: "${study_name}_chr${chr}.variant_id_map.tsv.gz", overwrite: true

    input:
    tuple val(study_name), val(chr), path(stage2_vcf), path(stage2_tbi), path(dbsnp_vcf), path(dbsnp_tbi)

    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}.pgen"), path("${study_name}_chr${chr}.pvar"), path("${study_name}_chr${chr}.psam"), path("${study_name}_chr${chr}.variant_qc.tsv"), path("${study_name}_chr${chr}.variant_id_map.tsv.gz")

    script:
    """
    INPUT_VARIANTS=\$(\$BCFTOOLS_BIN view -H "${stage2_vcf}" | wc -l | tr -d ' ')

    \$BCFTOOLS_BIN annotate -a "${dbsnp_vcf}" -c ID -Oz -o ${study_name}_chr${chr}.annotated.vcf.gz "${stage2_vcf}"

    \$PYTHON3_BIN "${projectDir}/bin/normalize_variant_ids.py" \\
      --input ${study_name}_chr${chr}.annotated.vcf.gz \\
      --output ${study_name}_chr${chr}.normalized.vcf \\
      --mapping-output ${study_name}_chr${chr}.variant_id_map.tsv.gz

    bgzip -f ${study_name}_chr${chr}.normalized.vcf
    \$BCFTOOLS_BIN index -f -t ${study_name}_chr${chr}.normalized.vcf.gz

    \$BCFTOOLS_BIN view \\
      -i 'INFO/R2>=${params.min_r2} && INFO/MAF>=${params.maf}' \\
      ${study_name}_chr${chr}.normalized.vcf.gz \\
      -Oz -o ${study_name}_chr${chr}.filtered.vcf.gz

    \$BCFTOOLS_BIN index -f -t ${study_name}_chr${chr}.filtered.vcf.gz

    if [ "${chr}" = "X" ]; then
      \$PLINK2_BIN \\
        --vcf ${study_name}_chr${chr}.filtered.vcf.gz dosage=HDS \\
        --double-id \\
        --split-par b38 \\
        --lax-chrx-import \\
        --make-pgen \\
        --sort-vars \\
        --out ${study_name}_chr${chr}.import
    else
      \$PLINK2_BIN \\
        --vcf ${study_name}_chr${chr}.filtered.vcf.gz dosage=HDS \\
        --double-id \\
        --make-pgen \\
        --sort-vars \\
        --out ${study_name}_chr${chr}.import
    fi

    POST_R2_MAF_VARIANTS=\$(grep -vc '^#' ${study_name}_chr${chr}.import.pvar)

    HWE_APPLIED=0
    if [ "${params.run_hwe}" = "true" ] && [ "${chr}" != "X" ]; then
      \$PLINK2_BIN \\
        --pfile ${study_name}_chr${chr}.import \\
        --hwe ${params.hwe_p} ${params.hwe_k} midp keep-fewhet \\
        --make-pgen \\
        --out ${study_name}_chr${chr}
      HWE_APPLIED=1
    else
      \$PLINK2_BIN \\
        --pfile ${study_name}_chr${chr}.import \\
        --make-pgen \\
        --out ${study_name}_chr${chr}
    fi

    FINAL_VARIANTS=\$(grep -vc '^#' ${study_name}_chr${chr}.pvar)

    read -r RSID_COUNT FALLBACK_COUNT DUPLICATE_RSID_FALLBACK_COUNT <<EOF_COUNTS
\$(gzip -cd ${study_name}_chr${chr}.variant_id_map.tsv.gz | awk -F'\\t' '
  NR == 1 { next }
  \$7 == "rsid" { rsid++ }
  \$7 == "fallback" { fallback++ }
  \$7 == "duplicate_rsid_fallback" { dup++; fallback++ }
  END { printf "%d %d %d\\n", rsid + 0, fallback + 0, dup + 0 }
')
EOF_COUNTS

    {
      printf "study\\tchr\\tinput_variants\\trsid_variants\\tfallback_variants\\tduplicate_rsid_fallbacks\\tpost_r2_maf_variants\\tfinal_variants\\thwe_applied\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${study_name}" "${chr}" "\${INPUT_VARIANTS}" "\${RSID_COUNT}" "\${FALLBACK_COUNT}" "\${DUPLICATE_RSID_FALLBACK_COUNT}" "\${POST_R2_MAF_VARIANTS}" "\${FINAL_VARIANTS}" "\${HWE_APPLIED}"
    } > ${study_name}_chr${chr}.variant_qc.tsv
    """
}
