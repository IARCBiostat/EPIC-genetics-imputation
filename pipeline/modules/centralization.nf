process CENTRALIZATION {
    label 'mid_mem'
    publishDir "${params.outdir}/centralized", mode: 'copy'

    input:
    path study_bfiles // List of [bed, bim, fam] for each study

    output:
    path "epic_centralized.bed"
    path "epic_centralized.bim"
    path "epic_centralized.fam"

    script:
    """
    # Logic to merge multiple studies and perform final cross-batch checks
    # For now, just a placeholder
    echo "Merging studies..."
    touch epic_centralized.bed epic_centralized.bim epic_centralized.fam
    """
}
