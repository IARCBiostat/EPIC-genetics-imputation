process NORMALIZATION {
    tag "${study_id}"
    label 'mid_mem'
    publishDir "${params.outdir}/${study_id}/normalization", mode: 'copy'

    input:
    tuple val(study_id), path(bed), path(bim), path(fam), path(manifest), path(id_linkage_file)

    output:
    tuple val(study_id), path("${study_id}_normalized.bed"), path("${study_id}_normalized.bim"), path("${study_id}_normalized.fam"), emit: normalized_data
    path "${study_id}_qc_report.txt", emit: report

    script:
    """
    # 1. Name Linkage (Bim vs Manifest)
    python2.7 ${baseDir}/bin/name_linkage_Bim_Manifest.py --file ${study_id} --man ${manifest}
    
    # 2. Initial QC
    python2.7 ${baseDir}/bin/preProcessing.py --file ${study_id} --man ${manifest}
    
    # 3. ID Standardization
    python2.7 ${baseDir}/bin/link_ID_standard.py --file ${study_id} --id ${id_linkage_file}
    
    # 4. Coordinate Completion
    python2.7 ${baseDir}/bin/search_chr_standard.py --file ${study_id} --man ${manifest}
    python2.7 ${baseDir}/bin/search_pos_standard.py --file ${study_id} --man ${manifest}
    python2.7 ${baseDir}/bin/search_alleles_standard.py --file ${study_id} --man ${manifest}
    
    # Final conversion to standard PLINK format for next steps
    # (Mocking final step for now)
    cp ${bed} ${study_id}_normalized.bed
    cp ${bim} ${study_id}_normalized.bim
    cp ${fam} ${study_id}_normalized.fam
    echo "Normalization complete." > ${study_id}_qc_report.txt
    """
}
