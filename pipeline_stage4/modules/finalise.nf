process FINALISE_STUDY {
    tag "${study_name}"
    publishDir "${params.dest_root}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "${study_name}/${filename}"
    }

    input:
    tuple val(study_name),
          path(pgens),
          path(pvars),
          path(psam),
          path(related_exclude),
          path(ancestry_exclude),
          path(hwe_exclude),
          path(stage3_html),
          path(master_html),
          path(het),
          path(kin0),
          path(eigenvec),
          path(eigenval)

    output:
    tuple val(study_name),
          path("${study_name}.tar.gz"),
          path("report-stage2.html"),
          path("report-stage3.html"),
          path("report-master.html"),
          path("review")

    script:
    def pgen_list = pgens.collect { it.toString() }.join(' ')
    def pvar_list = pvars.collect { it.toString() }.join(' ')
    """
    # stageInMode='copy' means pgen/pvar files are already local copies in the work dir.
    # Move them into the archive directory (no second copy needed).
    mkdir -p "${study_name}"

    for f in \$(echo "${pgen_list} ${pvar_list}" | tr ' ' '\\n' | sort -V); do
      mv "\$f" "${study_name}/"
    done

    # Single psam for all chromosomes
    cp "${psam}" "${study_name}/${study_name}.psam"

    _pigz=\$(command -v pigz 2>/dev/null || true)
    if [ -n "\$_pigz" ]; then
        tar -I "pigz -p ${task.cpus}" -cf "${study_name}.tar.gz" "${study_name}"
    else
        tar -czf "${study_name}.tar.gz" "${study_name}"
    fi

    # Review files — published alongside the archive, outside it
    mkdir -p review
    cp "${related_exclude}" review/related.exclude
    cp "${ancestry_exclude}" review/ancestry.exclude
    cp "${hwe_exclude}" review/hwe.exclude
    cp "${het}"     review/het.het
    cp "${kin0}"    review/king.kin0
    cp "${eigenvec}" review/pca.eigenvec
    cp "${eigenval}" review/pca.eigenval

    # HTML reports — published alongside the archive
    _s2html="${params.stage2_root}/${study_name}/stage2/report/report-stage2.html"
    if [ -f "\${_s2html}" ]; then
        cp "\${_s2html}" report-stage2.html
    else
        printf '<html><body><p>Stage 2 report not available for %s.</p></body></html>\n' "${study_name}" > report-stage2.html
    fi
    cp "${master_html}" report-master.html
    """
}
