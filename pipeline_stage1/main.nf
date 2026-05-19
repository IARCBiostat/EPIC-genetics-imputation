#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * EPIC Genetics Stage 1
 *
 * The study-specific Python scripts remain the authoritative implementation of
 * the filtering and harmonisation logic. Nextflow owns study-level orchestration
 * and report generation so each requested study can be scheduled independently.
 */

def studyFromScript(scriptPath) {
    def scriptFile = new File(scriptPath.toString())
    def matcher = scriptFile.text =~ /(?m)^STUDY_ID\s*=\s*['"]([^'"]+)['"]/
    if (!matcher.find()) {
        throw new IllegalArgumentException("Could not find STUDY_ID in ${scriptPath}")
    }
    return matcher.group(1)
}

def includedStudies() {
    params.study == 'all' ? [] : params.study.split(',').collect { it.trim() }.findAll { it }
}

process STAGE1_RAW_COPY {
    tag "${study}"

    input:
    tuple val(study), val(script_path)

    output:
    tuple val(study), val(script_path), path("${study}.raw_copy.done")

    script:
    """
    printf "study\\tstep\\n%s\\traw_copy\\n" "${study}" > ${study}.raw_copy.done
    """
}

process STAGE1_INITIAL_REVIEW {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.initial_review.done")

    script:
    """
    printf "study\\tstep\\n%s\\tinitial_review\\n" "${study}" > ${study}.initial_review.done
    """
}

process STAGE1_OPTIONAL_PRE_STEP {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.optional_pre_step.done")

    script:
    """
    printf "study\\tstep\\n%s\\toptional_pre_step\\n" "${study}" > ${study}.optional_pre_step.done
    """
}

process STAGE1_ID_LINKAGE {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.id_linkage.done")

    script:
    """
    printf "study\\tstep\\n%s\\tid_linkage\\n" "${study}" > ${study}.id_linkage.done
    """
}

process STAGE1_PART1_FILTERING {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.part1.done")

    script:
    """
    printf "study\\tstep\\n%s\\tpart1_filtering\\n" "${study}" > ${study}.part1.done
    """
}

process STAGE1_COMPLETION {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.completion.done")

    script:
    """
    printf "study\\tstep\\n%s\\tcompletion\\n" "${study}" > ${study}.completion.done
    """
}

process STAGE1_POST_COMPLETION_REVIEW {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.post_completion_review.done")

    script:
    """
    printf "study\\tstep\\n%s\\tpost_completion_review\\n" "${study}" > ${study}.post_completion_review.done
    """
}

process STAGE1_PART2_FILTERING {
    tag "${study}"

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), val(script_path), path("${study}.part2.done")

    script:
    """
    printf "study\\tstep\\n%s\\tpart2_filtering\\n" "${study}" > ${study}.part2.done
    """
}

process STAGE1_LIFTOVER_HANDOFF {
    tag "${study}"
    cpus 2
    memory '32 GB'

    input:
    tuple val(study), val(script_path), path(previous_step)

    output:
    tuple val(study), path("${study}.stage1.done")

    script:
    """
    ${params.python3_bin} "${script_path}" \\
      --data-root "${params.data_root}" \\
      --work-root "${params.work_root}" \\
      --outdir "${params.outdir}" \\
      --plink "${params.plink_bin}" \\
      --python2 "${params.python2_bin}"

    test -f "${params.outdir}/${study}/stage1/${study}.bed"
    test -f "${params.outdir}/${study}/stage1/${study}.bim"
    test -f "${params.outdir}/${study}/stage1/${study}.fam"
    printf "study\\tstep\\n%s\\tliftover_handoff\\n" "${study}" > ${study}.stage1.done
    """
}

process STAGE1_QC {
    tag "${study}"
    cpus 4
    memory '32 GB'

    input:
    tuple val(study), path(previous_step)

    output:
    tuple val(study), path("${study}.stage1_qc.done")

    script:
    def threads = task.cpus
    def backoff_secs = (task.attempt - 1) * 30
    """
    [ ${backoff_secs} -gt 0 ] && sleep ${backoff_secs}

    INPUT="${params.outdir}/${study}/stage1/${study}"
    QC_DIR="${params.outdir}/${study}/stage1/qc"
    mkdir -p "\${QC_DIR}"

    N_INIT_SAMPLES=\$(wc -l < "\${INPUT}.fam" | tr -d ' ')
    N_INIT_VARS=\$(wc -l < "\${INPUT}.bim" | tr -d ' ')

    # --- Sex check (requires >= ${params.sex_check_min_variants} chrX variants) ---
    CHR_X_VARS=\$(awk '\$1=="X" || \$1=="23" {count++} END {print count+0}' "\${INPUT}.bim")
    if [ "\${CHR_X_VARS}" -ge ${params.sex_check_min_variants} ]; then
        if ${params.plink_bin} --bfile "\${INPUT}" \\
                --check-sex ${params.sex_check_f_female} ${params.sex_check_f_male} \\
                --out sex_check; then
            awk 'NR>1 && \$5=="PROBLEM" && (\$4==1 || \$4==2) {print \$1, \$2}' \\
                sex_check.sexcheck > sex_mismatch.id
        else
            : > sex_mismatch.id
        fi
    else
        : > sex_mismatch.id
    fi
    N_SEX_MISMATCH=\$(wc -l < sex_mismatch.id | tr -d ' ')

    # --- Heterozygosity outliers (LD-pruned autosomes) ---
    ${params.plink2_bin} --bfile "\${INPUT}" \\
        --chr 1-22 \\
        --indep-pairwise 1500 150 0.2 \\
        --threads ${threads} \\
        --out pruned

    ${params.plink_bin} --bfile "\${INPUT}" \\
        --extract pruned.prune.in \\
        --het \\
        --out het_check

    awk 'NR>1 {sum+=\$6; sumsq+=\$6^2; n++} END {
        if(n>1){mean=sum/n; v=(sumsq-n*mean^2)/(n-1); if(v<0)v=0; print mean, sqrt(v)}
        else{print 0,0}
    }' het_check.het > het_stats.txt
    read HET_MEAN HET_SD < het_stats.txt
    awk -v mean="\${HET_MEAN}" -v sd="\${HET_SD}" -v thr="${params.het_sd_threshold}" \\
        'NR>1 && (\$6 < mean-thr*sd || \$6 > mean+thr*sd) {print \$1, \$2}' \\
        het_check.het > het_outliers.id
    N_HET_OUTLIERS=\$(wc -l < het_outliers.id | tr -d ' ')

    # --- KING duplicate removal (cutoff ${params.king_duplicate_cutoff} = MZ twin / duplicate threshold) ---
    ${params.plink2_bin} --bfile "\${INPUT}" \\
        --extract pruned.prune.in \\
        --make-king-table \\
        --threads ${threads} \\
        --out king

    [ -f king.kin0 ] || : > king.kin0
    if [ -s king.kin0 ]; then
        ${params.plink2_bin} --king-cutoff-table king.kin0 ${params.king_duplicate_cutoff} --out king_dupes
    else
        : > king_dupes.king.cutoff.out.id
    fi
    [ -f king_dupes.king.cutoff.out.id ] || : > king_dupes.king.cutoff.out.id
    N_DUPLICATES=\$(wc -l < king_dupes.king.cutoff.out.id | tr -d ' ')

    # --- Merge exclusion lists ---
    cat sex_mismatch.id het_outliers.id king_dupes.king.cutoff.out.id | sort -u > samples_to_remove.id
    N_REMOVED_SAMPLES=\$(wc -l < samples_to_remove.id | tr -d ' ')

    # --- Apply sample exclusions, HWE (autosomes), and pre-phasing MAF filter ---
    REMOVE_OPT=""
    [ -s samples_to_remove.id ] && REMOVE_OPT="--remove samples_to_remove.id"

    ${params.plink2_bin} --bfile "\${INPUT}" \\
        \${REMOVE_OPT} \\
        --hwe ${params.hwe_p} midp \\
        --maf ${params.maf_phasing} \\
        --threads ${threads} \\
        --make-bed \\
        --out ${study}

    N_FINAL_SAMPLES=\$(wc -l < ${study}.fam | tr -d ' ')
    N_FINAL_VARS=\$(wc -l < ${study}.bim | tr -d ' ')

    # --- Write QC report ---
    {
        printf "study\\tn_initial_samples\\tn_sex_mismatch\\tn_het_outliers\\tn_duplicates\\tn_removed_samples\\tn_final_samples\\tn_initial_variants\\tn_final_variants\\n"
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "${study}" "\${N_INIT_SAMPLES}" "\${N_SEX_MISMATCH}" "\${N_HET_OUTLIERS}" \\
            "\${N_DUPLICATES}" "\${N_REMOVED_SAMPLES}" "\${N_FINAL_SAMPLES}" \\
            "\${N_INIT_VARS}" "\${N_FINAL_VARS}"
    } > "\${QC_DIR}/${study}.stage1_qc.tsv"

    # Overwrite stage1 BED/BIM/FAM with QC-filtered versions
    cp ${study}.bed "\${INPUT}.bed"
    cp ${study}.bim "\${INPUT}.bim"
    cp ${study}.fam "\${INPUT}.fam"

    printf "study\\tstep\\n%s\\tstage1_qc\\n" "${study}" > ${study}.stage1_qc.done
    """
}

process STAGE1_REPORTS {
    tag "${study}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(study), path(previous_step)

    output:
    tuple val(study), path("${study}.reports.done")

    script:
    """
    ${params.python3_bin} "${projectDir}/bin/stage1_reports.py" \\
      --analysis-root "${params.outdir}" \\
      --studies "${study}"
    printf "study\\tstep\\n%s\\treports\\n" "${study}" > ${study}.reports.done
    """
}

workflow {
    def included = includedStudies()

    ch_scripts = Channel
        .fromPath("${projectDir}/scripts/process_*.py")
        .map { script_path -> tuple(studyFromScript(script_path), script_path.toString()) }
        .filter { study, script_path -> params.study == 'all' || included.contains(study) }

    ch_raw = STAGE1_RAW_COPY(ch_scripts)
    ch_initial = STAGE1_INITIAL_REVIEW(ch_raw)
    ch_pre = STAGE1_OPTIONAL_PRE_STEP(ch_initial)
    ch_linked = STAGE1_ID_LINKAGE(ch_pre)
    ch_part1 = STAGE1_PART1_FILTERING(ch_linked)
    ch_completed = STAGE1_COMPLETION(ch_part1)
    ch_reviewed = STAGE1_POST_COMPLETION_REVIEW(ch_completed)
    ch_part2 = STAGE1_PART2_FILTERING(ch_reviewed)
    ch_handoff = STAGE1_LIFTOVER_HANDOFF(ch_part2)
    ch_qc = STAGE1_QC(ch_handoff)
    STAGE1_REPORTS(ch_qc)
}
