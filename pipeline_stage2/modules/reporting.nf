process REPORTING {
    tag "${study_name}"

    input:
    tuple val(study_name), val(reporting_inputs)

    output:
    path("reporting.done")

    script:
    def flattened_inputs = reporting_inputs.flatten()
    def input_count = flattened_inputs.size()
    """
    echo "Reporting trigger for study: ${study_name}"
    echo "Input file count: ${input_count}"

    if ! \$PYTHON3_BIN -c "import matplotlib, pandas" >/dev/null 2>&1; then
        _py=\$(find "\${NXF_CONDA_CACHE}" -maxdepth 4 -name "python3" 2>/dev/null | while IFS= read -r _p; do
            [ -x "\${_p}" ] || continue
            "\${_p}" -c "import matplotlib, pandas" 2>/dev/null && printf '%s\n' "\${_p}" && break
        done)
        [ -n "\$_py" ] && PYTHON3_BIN="\$_py"
        unset _py
    fi

    \$PYTHON3_BIN ${projectDir}/bin/run_stage2_reports.py \\
      --analysis-root "${params.outdir}" \\
      --stage1-root "${params.stage1_root}" \\
      --summary-output "${params.outdir}/stage2-summary.md" \\
      --studies "${study_name}" \\
      --chromosomes "${params.chromosomes}" \\
      --reference-dir "${params.ref_1000g_dir}" \\
      --min-r2-threshold "${params.summary_min_r2}" \\
      --high-quality-threshold "${params.summary_high_quality}"

    printf "step\\tstatus\\nreporting\\tcomplete\\n" > reporting.done
    """

    stub:
    """
    printf "step\\tstatus\\nreporting\\tstub\\n" > reporting.done
    """
}
