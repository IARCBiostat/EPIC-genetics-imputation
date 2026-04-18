// pipeline/modules/phasing.nf

process PHASE_AUTOSOMES {
    tag "${study_name}/chr${chr}"
    cpus { params.eagle_threads }
    memory '32 GB'

    input:
    tuple val(study_name), val(chr), path(target_vcf), path(target_vcf_index), path(ref_bcf), path(ref_bcf_index)
    
    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}_GxS.phased.vcf.gz"), path("${study_name}_chr${chr}_GxS.phased.vcf.gz.csi")

    script:
    """
    EAGLE_CMD="\$(command -v "\${EAGLE_BIN}" 2>/dev/null || command -v eagle 2>/dev/null || command -v Eagle 2>/dev/null || true)"
    if [ -z "\${EAGLE_CMD}" ]; then
        echo "Eagle binary not found in process environment" >&2
        exit 1
    fi

    "\${EAGLE_CMD}" \\
        --vcfRef ${ref_bcf} \\
        --vcfTarget ${target_vcf} \\
        --allowRefAltSwap \\
        --geneticMapFile=${params.genetic_map_file} \\
        --numThreads ${params.eagle_threads} \\
        --outPrefix ${study_name}_chr${chr}_GxS.phased

    \$BCFTOOLS_BIN index ${study_name}_chr${chr}_GxS.phased.vcf.gz
    """
}
