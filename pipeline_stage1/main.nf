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
    STAGE1_REPORTS(ch_handoff)
}
