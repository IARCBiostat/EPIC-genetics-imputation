process PREP_TARGET_VCF {
    tag "${study_name}/chr${chr}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename ->
        def output_name = filename.tokenize('/').last()
        def study = output_name.replaceFirst(/_chr.*$/, '')
        return "${study}/stage2/target/${output_name}"
    }

    input:
    tuple val(study_name), val(chr), val(base), path(bed), path(bim), path(fam)

    output:
    tuple val(chr), val(study_name), path("target/${study_name}_chr${chr}.vcf.gz"), path("target/${study_name}_chr${chr}.vcf.gz.tbi"), optional: true

    script:
    """
    mkdir -p target

    cat > chr_name_conv.txt <<'EOF'
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
18	chr18
19	chr19
20	chr20
21	chr21
22	chr22
X	chrX
Y	chrY
23	chrX
24	chrY
25	chrXY
26	chrMT
MT	chrMT
M	chrMT
EOF

    if ! awk -v chr="${chr}" '
        \$1 == chr || (chr == "X" && (\$1 == "23" || \$1 == "chrX")) {
            found = 1
            exit 0
        }
        END {
            exit(found ? 0 : 1)
        }
    ' ${bim} ; then
        echo "Skipping chr${chr}: no variants present in the stage-1 BIM"
        exit 0
    fi

    \$PLINK_BIN --bfile ${bed.baseName} \\
          --chr ${chr} \\
          --recode vcf bgz \\
          --out tmp_${study_name}_chr${chr}

    \$BCFTOOLS_BIN annotate --rename-chrs chr_name_conv.txt tmp_${study_name}_chr${chr}.vcf.gz -Ou | \\
        \$BCFTOOLS_BIN annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou | \\
        \$BCFTOOLS_BIN norm -m -any -Ou | \\
        \$BCFTOOLS_BIN +fill-tags -Oz -o target/${study_name}_chr${chr}.vcf.gz -- -t AC,AN

    \$BCFTOOLS_BIN index -f -t target/${study_name}_chr${chr}.vcf.gz

    rm -f tmp_${study_name}_chr${chr}.log tmp_${study_name}_chr${chr}.nosex
    rm -f tmp_${study_name}_chr${chr}.vcf.gz tmp_${study_name}_chr${chr}.vcf.gz.tbi
    """
}
