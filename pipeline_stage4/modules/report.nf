process MASTER_REPORT {
    tag "${study_name}"
    publishDir "${params.report_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(study_name), path(stage3_html)

    output:
    tuple val(study_name), path("${study_name}.master-report.html")

    script:
    """
    \$PYTHON3_BIN "${params.project_root}/src/misc/master_report.py" \\
      --analysis-root ${params.analysis_root} \\
      --report-dir . \\
      --studies ${study_name}
    """
}
