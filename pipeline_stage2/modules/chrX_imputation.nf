// pipeline/modules/chrX_imputation.nf

process IMPUTE_CHRX {
    tag "${study_name}/chrX"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        def study = filename.replaceFirst(/_chr.*$/, '')
        return "${study}/stage2/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(target_vcf), path(target_tbi), path(msavs), path(bcfs)

    output:
    tuple val(study_name), path("${study_name}_chrX_GxS.imputed.vcf.gz"), path("${study_name}_chrX_GxS.imputed.vcf.gz.tbi")

    script:
    def threads = task.cpus
    """
    set -euo pipefail

    prepare_minimac_target() {
        local block="\$1"
        local input="phased.\${block}.vcf.gz"
        local output="phased.\${block}.minimac.vcf.gz"

        if [ "${params.chrx_force_diploid}" = "true" ]; then
            if ! \$BCFTOOLS_BIN +fixploidy "\${input}" -Oz -o "\${output}" -- -f 2; then
                echo "Could not normalize chrX \${block} ploidy for Minimac4" >&2
                return 1
            fi
            \$BCFTOOLS_BIN index -f -t "\${output}"
        else
            cp "\${input}" "\${output}"
            cp "\${input}.tbi" "\${output}.tbi"
        fi
    }

    echo -e "chrX\\t2531036\\t2531036\\tHGSV_248479\\nchrX\\t2777488\\t2777488\\tHGSV_249395" > exclude_ids.bed
    declare -a imputed_blocks=()

    # ---------------------------------------------------------
    # 1. PAR1 Block (chrX:1-2781513)
    # ---------------------------------------------------------
    \$BCFTOOLS_BIN query -f '%CHROM\\t%POS\\t%POS\\n' 1kGP_panel_chrX.PAR1.bcf > ref.PAR1.positions.tsv
    \$BCFTOOLS_BIN view -r chrX:1-2781513 ${target_vcf} -Oz -o target.PAR1.raw.vcf.gz
    \$BCFTOOLS_BIN index -t target.PAR1.raw.vcf.gz
    \$BCFTOOLS_BIN view -R ref.PAR1.positions.tsv target.PAR1.raw.vcf.gz -Oz -o target.PAR1.vcf.gz
    \$BCFTOOLS_BIN index -t target.PAR1.vcf.gz

    if [ "\$(\$BCFTOOLS_BIN view -H target.PAR1.vcf.gz | wc -l | tr -d ' ')" -gt 0 ]; then
        if \$SHAPEIT5_COMMON_BIN \\
            --input target.PAR1.vcf.gz \\
            --reference 1kGP_panel_chrX.PAR1.bcf \\
            --map ${params.shapeit5_map_dir}/chrX.b38.gmap.gz \\
            --region chrX:1-2781513 \\
            --thread ${threads} \\
            --output phased.PAR1.vcf.gz; then
            \$BCFTOOLS_BIN index -t phased.PAR1.vcf.gz

            if prepare_minimac_target PAR1 && \\
               \$MINIMAC4_BIN 1kGP_panel_chrX.PAR1.msav phased.PAR1.minimac.vcf.gz \\
                -o imputed.PAR1.vcf.gz --threads ${params.minimac_threads} -b ${params.minimac_batch_size} --min-r2 ${params.min_r2} -O vcf.gz; then
                \$BCFTOOLS_BIN index -t imputed.PAR1.vcf.gz
                imputed_blocks+=(imputed.PAR1.vcf.gz)
            else
                echo "Skipping chrX PAR1: Minimac4 could not impute this block"
                rm -f imputed.PAR1.vcf.gz imputed.PAR1.vcf.gz.tbi
            fi
        else
            echo "Skipping chrX PAR1: SHAPEIT5 failed on this block"
            rm -f phased.PAR1.vcf.gz phased.PAR1.vcf.gz.tbi
        fi
    else
        echo "Skipping chrX PAR1: no overlapping target/reference variants after filtering"
    fi

    # ---------------------------------------------------------
    # 2. nonPAR Block (chrX:2781514-155700882)
    # ---------------------------------------------------------
    \$BCFTOOLS_BIN query -f '%CHROM\\t%POS\\t%POS\\n' 1kGP_panel_chrX.nonPAR.bcf > ref.nonPAR.positions.tsv
    \$BCFTOOLS_BIN view -r chrX:2781514-155700882 -T ^exclude_ids.bed -i 'POS>=2781514' ${target_vcf} -Oz -o target.nonPAR.raw.vcf.gz
    \$BCFTOOLS_BIN index -t target.nonPAR.raw.vcf.gz
    \$BCFTOOLS_BIN view -R ref.nonPAR.positions.tsv target.nonPAR.raw.vcf.gz -Oz -o target.nonPAR.vcf.gz
    \$BCFTOOLS_BIN index -t target.nonPAR.vcf.gz

    if [ "\$(\$BCFTOOLS_BIN view -H target.nonPAR.vcf.gz | wc -l | tr -d ' ')" -gt 0 ]; then
        if \$SHAPEIT5_COMMON_BIN \\
            --input target.nonPAR.vcf.gz \\
            --reference 1kGP_panel_chrX.nonPAR.bcf \\
            --map ${params.shapeit5_map_dir}/chrX.b38.gmap.gz \\
            --region chrX:2781514-155700882 \\
            --thread ${threads} \\
            --output phased.nonPAR.vcf.gz; then
            \$BCFTOOLS_BIN index -t phased.nonPAR.vcf.gz

            if prepare_minimac_target nonPAR && \\
               \$MINIMAC4_BIN 1kGP_panel_chrX.nonPAR.msav phased.nonPAR.minimac.vcf.gz \\
                -o imputed.nonPAR.vcf.gz --threads ${params.minimac_threads} -b ${params.minimac_batch_size} --min-ratio ${params.chrx_min_ratio} --min-r2 ${params.min_r2} -O vcf.gz; then
                \$BCFTOOLS_BIN index -t imputed.nonPAR.vcf.gz
                imputed_blocks+=(imputed.nonPAR.vcf.gz)
            else
                echo "Skipping chrX nonPAR: Minimac4 could not impute this block"
                rm -f imputed.nonPAR.vcf.gz imputed.nonPAR.vcf.gz.tbi
            fi
        else
            echo "Skipping chrX nonPAR: SHAPEIT5 failed on this block"
            rm -f phased.nonPAR.vcf.gz phased.nonPAR.vcf.gz.tbi
        fi
    else
        echo "Skipping chrX nonPAR: no overlapping target/reference variants after filtering"
    fi

    # ---------------------------------------------------------
    # 3. PAR2 Block (chrX:155700883-)
    # ---------------------------------------------------------
    \$BCFTOOLS_BIN query -f '%CHROM\\t%POS\\t%POS\\n' 1kGP_panel_chrX.PAR2.bcf > ref.PAR2.positions.tsv
    \$BCFTOOLS_BIN view -r chrX:155700883- -i 'POS>=155700883' ${target_vcf} -Oz -o target.PAR2.raw.vcf.gz
    \$BCFTOOLS_BIN index -t target.PAR2.raw.vcf.gz
    \$BCFTOOLS_BIN view -R ref.PAR2.positions.tsv target.PAR2.raw.vcf.gz -Oz -o target.PAR2.vcf.gz
    \$BCFTOOLS_BIN index -t target.PAR2.vcf.gz

    if [ "\$(\$BCFTOOLS_BIN view -H target.PAR2.vcf.gz | wc -l | tr -d ' ')" -gt 0 ]; then
        if \$SHAPEIT5_COMMON_BIN \\
            --input target.PAR2.vcf.gz \\
            --reference 1kGP_panel_chrX.PAR2.bcf \\
            --map ${params.shapeit5_map_dir}/chrX.b38.gmap.gz \\
            --region chrX:155700883- \\
            --thread ${threads} \\
            --output phased.PAR2.vcf.gz; then
            \$BCFTOOLS_BIN index -t phased.PAR2.vcf.gz

            if prepare_minimac_target PAR2 && \\
               \$MINIMAC4_BIN 1kGP_panel_chrX.PAR2.msav phased.PAR2.minimac.vcf.gz \\
                -o imputed.PAR2.vcf.gz --threads ${params.minimac_threads} -b ${params.minimac_batch_size} --min-r2 ${params.min_r2} -O vcf.gz; then
                \$BCFTOOLS_BIN index -t imputed.PAR2.vcf.gz
                imputed_blocks+=(imputed.PAR2.vcf.gz)
            else
                echo "Skipping chrX PAR2: Minimac4 could not impute this block"
                rm -f imputed.PAR2.vcf.gz imputed.PAR2.vcf.gz.tbi
            fi
        else
            echo "Skipping chrX PAR2: SHAPEIT5 failed on this block"
            rm -f phased.PAR2.vcf.gz phased.PAR2.vcf.gz.tbi
        fi
    else
        echo "Skipping chrX PAR2: no overlapping target/reference variants after filtering"
    fi

    # ---------------------------------------------------------
    # 4. Concatenate and Annotate
    # ---------------------------------------------------------
    if [ "\${#imputed_blocks[@]}" -eq 0 ]; then
        \$BCFTOOLS_BIN view -h ${target_vcf} > empty.chrX.vcf
        \$BCFTOOLS_BIN view -Oz -o ${study_name}_chrX_GxS.imputed.vcf.gz empty.chrX.vcf
        \$BCFTOOLS_BIN index -t ${study_name}_chrX_GxS.imputed.vcf.gz
    elif [ "\${#imputed_blocks[@]}" -eq 1 ]; then
        cp "\${imputed_blocks[0]}" imputed.chrX.all.vcf.gz
        cp "\${imputed_blocks[0]}.tbi" imputed.chrX.all.vcf.gz.tbi
        \$BCFTOOLS_BIN annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' \\
            -Oz -o ${study_name}_chrX_GxS.imputed.vcf.gz \\
            imputed.chrX.all.vcf.gz
        \$BCFTOOLS_BIN index -t ${study_name}_chrX_GxS.imputed.vcf.gz
    else
        \$BCFTOOLS_BIN concat \\
            --threads 4 \\
            -Oz -o imputed.chrX.all.vcf.gz \\
            "\${imputed_blocks[@]}"
        \$BCFTOOLS_BIN index -t imputed.chrX.all.vcf.gz

        \$BCFTOOLS_BIN annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' \\
            -Oz -o ${study_name}_chrX_GxS.imputed.vcf.gz \\
            imputed.chrX.all.vcf.gz
        \$BCFTOOLS_BIN index -t ${study_name}_chrX_GxS.imputed.vcf.gz
    fi

    """
}
