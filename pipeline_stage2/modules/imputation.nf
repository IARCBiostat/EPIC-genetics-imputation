// pipeline/modules/imputation.nf

process IMPUTE_AUTOSOMES {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        def study = filename.replaceFirst(/_chr.*$/, '')
        return "${study}/stage2/${filename}"
    }

    input:
    tuple val(study_name), val(chr), path(phased_vcf), path(phased_vcf_index), path(msav_panel)
    
    output:
    tuple val(study_name), val(chr), path("${study_name}_chr${chr}_GxS.imputed.vcf.gz"), path("${study_name}_chr${chr}_GxS.imputed.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    output="${study_name}_chr${chr}_GxS.imputed.vcf.gz"

    run_minimac() {
        local stdout_log="\$1"
        local stderr_log="\$2"
        shift 2
        \$MINIMAC4_BIN \\
            ${msav_panel} \\
            ${phased_vcf} \\
            -o "\$output" \\
            --threads ${params.minimac_threads} \\
            -b ${params.minimac_batch_size} \\
            --min-r2 ${params.min_r2} \\
            -O vcf.gz \\
            "\$@" \\
            > "\$stdout_log" 2> "\$stderr_log"
    }

    set +e
    run_minimac minimac.full.stdout.log minimac.full.stderr.log
    minimac_status=\$?
    set -e

    if [ "\$minimac_status" -ne 0 ]; then
        echo "Minimac4 full-chromosome run failed for ${study_name}/chr${chr}; retrying over the target marker span." >&2
        cat minimac.full.stdout.log >&2 || true
        cat minimac.full.stderr.log >&2 || true
        rm -f "\$output" "\$output.tbi"

        read -r target_start target_end < <(\$BCFTOOLS_BIN query -f '%POS\\n' ${phased_vcf} | awk 'NR==1{first=\$1} {last=\$1} END{if(NR==0) exit 1; print first, last}')
        region="chr${chr}:\${target_start}-\${target_end}"
        help_text="\$(\$MINIMAC4_BIN --help 2>&1 || true)"

        set +e
        if printf '%s\\n' "\$help_text" | grep -q -- '--region'; then
            echo "Retrying Minimac4 with --region \$region" >&2
            run_minimac minimac.region.stdout.log minimac.region.stderr.log --region "\$region"
        elif printf '%s\\n' "\$help_text" | grep -q -- '--start' && printf '%s\\n' "\$help_text" | grep -q -- '--end'; then
            echo "Retrying Minimac4 with --chr ${chr} --start \$target_start --end \$target_end --window 0" >&2
            run_minimac minimac.region.stdout.log minimac.region.stderr.log --chr ${chr} --start "\$target_start" --end "\$target_end" --window 0
        elif printf '%s\\n' "\$help_text" | grep -q -- '--from' && printf '%s\\n' "\$help_text" | grep -q -- '--to'; then
            echo "Retrying Minimac4 with --chr ${chr} --from \$target_start --to \$target_end --window 0" >&2
            run_minimac minimac.region.stdout.log minimac.region.stderr.log --chr ${chr} --from "\$target_start" --to "\$target_end" --window 0
        else
            echo "Retrying Minimac4 with --region \$region; minimac4 --help did not expose a recognized region syntax." >&2
            run_minimac minimac.region.stdout.log minimac.region.stderr.log --region "\$region"
        fi
        minimac_status=\$?
        set -e

        if [ "\$minimac_status" -ne 0 ]; then
            cat minimac.region.stdout.log >&2 || true
            cat minimac.region.stderr.log >&2 || true
            exit "\$minimac_status"
        fi
    fi

    \$BCFTOOLS_BIN index -t "\$output"
    """
}
