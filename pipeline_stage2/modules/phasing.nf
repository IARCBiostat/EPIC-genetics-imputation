// pipeline/modules/phasing.nf

process PHASE_AUTOSOMES {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${filename.replaceFirst(/_chr.*$/, '')}/stage2/phasing/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(target_vcf), path(target_vcf_index), path(ref_bcf), path(ref_bcf_index)
    
    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}_GxS.phased.vcf.gz"), path("${study_name}_chr${chr}_GxS.phased.vcf.gz.csi")

    script:
    def threads = task.cpus
    """
    # SHAPEIT5 only recognises .bcf output; convert to vcf.gz afterwards
    \$SHAPEIT5_COMMON_BIN \\
        --input ${target_vcf} \\
        --reference ${ref_bcf} \\
        --map ${params.shapeit5_map_dir}/chr${chr}.b38.gmap.gz \\
        --region chr${chr} \\
        --thread ${threads} \\
        --output ${study_name}_chr${chr}_GxS.phased.bcf

    \$BCFTOOLS_BIN view -Oz -o ${study_name}_chr${chr}_GxS.phased.vcf.gz \\
        ${study_name}_chr${chr}_GxS.phased.bcf

    \$BCFTOOLS_BIN index ${study_name}_chr${chr}_GxS.phased.vcf.gz
    """
}
