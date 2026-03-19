process LIFTOVER {
    tag "${study_id}"
    label 'mid_mem'
    publishDir "${params.outdir}/${study_id}/liftover", mode: 'copy'

    input:
    tuple val(study_id), path(bed), path(bim), path(fam), path(chain_file)

    output:
    tuple val(study_id), path("${study_id}_hg38.bed"), path("${study_id}_hg38.bim"), path("${study_id}_hg38.fam"), emit: lifted_data

    script:
    """
    # Liftover logic using R and chain files
    # Mocking for now
    cp ${bed} ${study_id}_hg38.bed
    cp ${bim} ${study_id}_hg38.bim
    cp ${fam} ${study_id}_hg38.fam
    """
}
