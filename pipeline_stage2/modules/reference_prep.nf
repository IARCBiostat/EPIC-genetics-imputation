// pipeline/modules/reference_prep.nf

process PREP_REFERENCE {
    tag "chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        "cohort/stage2/reference/${filename}"
    }

    input:
    val chr

    output:
    tuple val(chr), path("*.msav"), path("*.bcf*")

    script:
    def vcf_prefix = chr == 'X' ? 
        "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2" : 
        "1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel"
    
    if ( chr != 'X' )
        """
        # 1. Convert to BCF, filter MAC < 10, normalize against fasta
        \$BCFTOOLS_BIN view --no-version -c 10 \\
            ${params.ref_1000g_dir}/${vcf_prefix}.vcf.gz | \\
        \$BCFTOOLS_BIN norm --no-version -Ou -m -any | \\
        \$BCFTOOLS_BIN norm --no-version -Ob \\
            -o 1kGP_panel_chr${chr}.bcf \\
            -d none -f ${params.fasta_ref}
            
        \$BCFTOOLS_BIN index -f 1kGP_panel_chr${chr}.bcf

        # 2. Compress to MSAV format for Minimac4
        \$MINIMAC4_BIN --compress-reference 1kGP_panel_chr${chr}.bcf > 1kGP_panel_chr${chr}.msav
        """
    else
        """
        # Minimac4 segfaulted on chrX PAR2 when fed SV records from the mixed
        # SNV/INDEL/SV panel, so chrX reference prep keeps only SNPs and indels.

        # 1. Prep PAR1 (1-2781513)
        \$BCFTOOLS_BIN view --no-version -c 10 -v snps,indels -r chrX:1-2781513 ${params.ref_1000g_dir}/${vcf_prefix}.vcf.gz | \\
        \$BCFTOOLS_BIN norm --no-version -Ou -m -any | \\
        \$BCFTOOLS_BIN norm --no-version -Ob -o 1kGP_panel_chrX.PAR1.bcf -d none -f ${params.fasta_ref}
        \$BCFTOOLS_BIN index -f 1kGP_panel_chrX.PAR1.bcf
        \$MINIMAC4_BIN --compress-reference 1kGP_panel_chrX.PAR1.bcf > 1kGP_panel_chrX.PAR1.msav

        # 2. Prep nonPAR (2781514-155700882)
        \$BCFTOOLS_BIN view --no-version -c 10 -v snps,indels -r chrX:2781514-155700882 ${params.ref_1000g_dir}/${vcf_prefix}.vcf.gz | \\
        \$BCFTOOLS_BIN norm --no-version -Ou -m -any | \\
        \$BCFTOOLS_BIN norm --no-version -Ob -o 1kGP_panel_chrX.nonPAR.bcf -d none -f ${params.fasta_ref}
        \$BCFTOOLS_BIN index -f 1kGP_panel_chrX.nonPAR.bcf
        \$MINIMAC4_BIN --compress-reference 1kGP_panel_chrX.nonPAR.bcf > 1kGP_panel_chrX.nonPAR.msav

        # 3. Prep PAR2 (155700883-)
        \$BCFTOOLS_BIN view --no-version -c 10 -v snps,indels -r chrX:155700883- ${params.ref_1000g_dir}/${vcf_prefix}.vcf.gz | \\
        \$BCFTOOLS_BIN norm --no-version -Ou -m -any | \\
        \$BCFTOOLS_BIN norm --no-version -Ob -o 1kGP_panel_chrX.PAR2.bcf -d none -f ${params.fasta_ref}
        \$BCFTOOLS_BIN index -f 1kGP_panel_chrX.PAR2.bcf
        \$MINIMAC4_BIN --compress-reference 1kGP_panel_chrX.PAR2.bcf > 1kGP_panel_chrX.PAR2.msav
        """
}
