process PREP_DBSNP_CHROM {
    tag "chr${chr}"

    input:
    tuple val(chr), path(dbsnp_vcf), path(dbsnp_tbi)

    output:
    tuple val(chr), path("dbsnp_chr${chr}.vcf.gz"), path("dbsnp_chr${chr}.vcf.gz.tbi")

    script:
    """
    case "${chr}" in
      1) NC_CONTIG="NC_000001.11" ;;
      2) NC_CONTIG="NC_000002.12" ;;
      3) NC_CONTIG="NC_000003.12" ;;
      4) NC_CONTIG="NC_000004.12" ;;
      5) NC_CONTIG="NC_000005.10" ;;
      6) NC_CONTIG="NC_000006.12" ;;
      7) NC_CONTIG="NC_000007.14" ;;
      8) NC_CONTIG="NC_000008.11" ;;
      9) NC_CONTIG="NC_000009.12" ;;
      10) NC_CONTIG="NC_000010.11" ;;
      11) NC_CONTIG="NC_000011.10" ;;
      12) NC_CONTIG="NC_000012.12" ;;
      13) NC_CONTIG="NC_000013.11" ;;
      14) NC_CONTIG="NC_000014.9" ;;
      15) NC_CONTIG="NC_000015.10" ;;
      16) NC_CONTIG="NC_000016.10" ;;
      17) NC_CONTIG="NC_000017.11" ;;
      18) NC_CONTIG="NC_000018.10" ;;
      19) NC_CONTIG="NC_000019.10" ;;
      20) NC_CONTIG="NC_000020.11" ;;
      21) NC_CONTIG="NC_000021.9" ;;
      22) NC_CONTIG="NC_000022.11" ;;
      X) NC_CONTIG="NC_000023.11" ;;
      *) echo "Unsupported chromosome: ${chr}" >&2; exit 1 ;;
    esac

    printf "%s\\tchr%s\\n" "\${NC_CONTIG}" "${chr}" > chr_name_map.tsv

    \$BCFTOOLS_BIN view -r "\${NC_CONTIG}" "${dbsnp_vcf}" -Ou | \\
      \$BCFTOOLS_BIN annotate --rename-chrs chr_name_map.tsv -Oz -o dbsnp_chr${chr}.vcf.gz

    \$BCFTOOLS_BIN index -f -t dbsnp_chr${chr}.vcf.gz
    """
}
