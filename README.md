# Human Genetics Data Imputation Pipeline
## Using 1000 Genomes Phase 3 Reference Panel
### HRC-Compatible QC Standards

**Author:** Pipeline Documentation  
**Date:** February 27, 2026  
**Purpose:** Complete guide for imputing human genetics data using publicly available 1000 Genomes reference panel with QC standards compatible with HRC (Haplotype Reference Consortium) imputation

**Note:** This pipeline uses 1000 Genomes Phase 3 (publicly available) but implements the same stringent QC parameters used in HRC-based imputation pipelines to ensure compatibility with collaborators using HRC reference data.

---

## Table of Contents

1. [Overview](#overview)
2. [Comparison with Collaborator Pipeline](#comparison-your-pipeline-vs-collaborators-hrc-pipeline)
3. [Software Requirements](#software-requirements)
4. [Step 1: Download Reference Panel](#step-1-download-1000-genomes-phase-3-reference-panel)
5. [Step 2: Prepare Your Data](#step-2-prepare-your-genotype-data)
6. [Step 3: Quality Control](#step-3-quality-control-critical)
7. [Step 4: Strand Alignment](#step-4-strand-alignment-and-allele-matching)
8. [Step 5: Pre-phasing](#step-5-pre-phasing-your-data)
9. [Step 6: Imputation](#step-6-imputation-with-minimac4-or-impute5)
10. [Step 7: Post-Imputation QC](#step-7-post-imputation-quality-control)
11. [Step 8: Merge Results](#step-8-merge-chromosomes-optional)
12. [Step 9: Liftover to GRCh38](#step-9-liftover-from-grch37-to-grch38)
13. [Quality Metrics](#key-quality-metrics-to-check)
14. [Computational Resources](#computational-resources)

---

## Overview

**Input:** Your genotyping data (SNP array data in PLINK or VCF format)  
**Output:** Imputed genotypes with quality scores  
**Reference:** 1000 Genomes Phase 3 (2,504 individuals from 26 populations)

**Important Note:** This pipeline is configured to match HRC (Haplotype Reference Consortium) imputation standards. While we use the publicly available 1000 Genomes Phase 3 reference panel instead of HRC (which requires data access application), the QC parameters and filtering thresholds match those used in HRC-based imputation pipelines.

**QC Parameters Match:**
- Pre-imputation: SNP call rate >98%, MAF >1%, HWE p>5×10⁻⁶, sample call rate >99%
- Post-imputation: R² >0.7, MAF >1%, HWE p>5×10⁻⁶, comprehensive sample QC
- These stringent thresholds ensure compatibility with collaborators using HRC reference panel

**Pipeline Summary:**
1. Download 1000 Genomes reference panel
2. Prepare and QC your genotype data (matching HRC standards)
3. Align strands and match alleles to reference
4. Phase your genotypes
5. Perform imputation
6. Filter by imputation quality (R² >0.7)
7. Comprehensive sample-level QC
8. Merge chromosomes
9. Liftover to GRCh38 (optional)

---

## Software Requirements

Install the following tools before beginning:

- **PLINK 1.9+** - QC and format conversion
  - Download: https://www.cog-genomics.org/plink/
- **bcftools** - VCF manipulation
  - Download: http://samtools.github.io/bcftools/
- **tabix** - Indexing VCF files
  - Download: http://www.htslib.org/
- **SHAPEIT4** or **Eagle2** - Phasing
  - SHAPEIT4: https://odelaneau.github.io/shapeit4/
  - Eagle2: https://alkesgroup.broadinstitute.org/Eagle/
- **Minimac4** or **IMPUTE5** - Imputation
  - Minimac4: https://genome.sph.umich.edu/wiki/Minimac4
  - IMPUTE5: https://jmarchini.org/impute-5/
- **conform-gt** - Strand alignment
  - Download: http://faculty.washington.edu/browning/conform-gt.html
- **Picard Tools** - Liftover (for GRCh37 to GRCh38 conversion)
  - Download: https://broadinstitute.github.io/picard/
- **samtools** - Reference genome indexing
  - Download: http://www.htslib.org/
- **Java 8+** - For conform-gt and Picard

---

## Step 1: Download 1000 Genomes Phase 3 Reference Panel

Create a directory structure and download the reference files:

```bash
# Create directory structure
mkdir -p 1000G_reference
cd 1000G_reference

# Download phased VCF files for all chromosomes
for chr in {1..22} X; do
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
done

# Download genetic maps for phasing
mkdir genetic_maps
cd genetic_maps

# Download genetic map files from:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/

# Download all chr*.b37.gmap.gz files
```

**Note:** Total download size is approximately 200-300 GB.

---

## Step 2: Prepare Your Genotype Data

Convert your data to VCF format and split by chromosome:

```bash
# If you have PLINK files (.bed/.bim/.fam), convert to VCF
plink --bfile your_data --recode vcf --out your_data

# Split by chromosome for parallel processing
for chr in {1..22}; do
    plink --bfile your_data --chr ${chr} --recode vcf --out your_data_chr${chr}
done
```

**Expected input formats:**
- PLINK binary (.bed/.bim/.fam)
- VCF (.vcf or .vcf.gz)

---

## Step 3: Quality Control (Critical!)

Perform stringent QC on your genotype data before imputation. These parameters match standard practices and align with HRC imputation requirements.

### Step 3A: SNP-Level QC (Pre-Imputation)

```bash
# First, perform SNP QC on the complete dataset (all chromosomes together)
plink --bfile your_data \
      --geno 0.02 \
      --maf 0.01 \
      --hwe 0.000005 \
      --make-bed \
      --out your_data_snp_qc

# Expected filters:
# --geno 0.02: Remove SNPs with genotyping rate < 98%
# --maf 0.01: Remove SNPs with MAF < 1%
# --hwe 0.000005: Remove SNPs violating HWE (p < 5×10⁻⁶) in controls
```

**Important Note on HWE filtering:**
If you have case-control data, apply HWE filter only to controls:
```bash
# Extract controls (assuming phenotype coded as 1=control, 2=case)
plink --bfile your_data \
      --filter-controls \
      --hwe 0.000005 \
      --write-snplist \
      --out control_hwe_pass

# Apply all filters
plink --bfile your_data \
      --extract control_hwe_pass.snplist \
      --geno 0.02 \
      --maf 0.01 \
      --make-bed \
      --out your_data_snp_qc
```

### Step 3B: Sample-Level QC (Pre-Imputation)

```bash
# 1. Remove samples with low genotyping rate (<99%)
plink --bfile your_data_snp_qc \
      --mind 0.01 \
      --make-bed \
      --out your_data_sample_qc1

# 2. Check for contaminated/related samples
plink --bfile your_data_sample_qc1 \
      --genome \
      --out your_data_relatedness

# Manually review .genome file for PI_HAT values
# Remove contaminated samples or one from each pair with PI_HAT > 0.185 (2nd degree relatives)
# Create a file "samples_to_remove.txt" with FID and IID of samples to exclude

plink --bfile your_data_sample_qc1 \
      --remove samples_to_remove.txt \
      --make-bed \
      --out your_data_sample_qc2

# 3. Check sex concordance
plink --bfile your_data_sample_qc2 \
      --check-sex \
      --out sex_check

# Remove samples where sex cannot be defined (STATUS=0 in sex_check.sexcheck)
# Also identify sex mismatches (PEDSEX != SNPSEX) for later removal
awk '$4 == 0 {print $1, $2}' sex_check.sexcheck > sex_problems.txt

plink --bfile your_data_sample_qc2 \
      --remove sex_problems.txt \
      --make-bed \
      --out your_data_qc_final
```

### Step 3C: Split by Chromosome for Imputation

```bash
# After completing all QC, split by chromosome
for chr in {1..22}; do
    plink --bfile your_data_qc_final \
          --chr ${chr} \
          --recode vcf \
          --out your_data_chr${chr}_qc
    
    # Compress and index
    bgzip your_data_chr${chr}_qc.vcf
    tabix -p vcf your_data_chr${chr}_qc.vcf.gz
done
```

### QC Summary Statistics

```bash
# Generate QC report
echo "=== Pre-Imputation QC Summary ===" > qc_report.txt
echo "" >> qc_report.txt

# Original sample/SNP counts
plink --bfile your_data --freq --out original_stats
orig_samples=$(wc -l < your_data.fam)
orig_snps=$(wc -l < your_data.bim)

# Final sample/SNP counts
final_samples=$(wc -l < your_data_qc_final.fam)
final_snps=$(wc -l < your_data_qc_final.bim)

echo "Original samples: ${orig_samples}" >> qc_report.txt
echo "Original SNPs: ${orig_snps}" >> qc_report.txt
echo "Final samples: ${final_samples}" >> qc_report.txt
echo "Final SNPs: ${final_snps}" >> qc_report.txt
echo "Samples removed: $((orig_samples - final_samples))" >> qc_report.txt
echo "SNPs removed: $((orig_snps - final_snps))" >> qc_report.txt

cat qc_report.txt
```

**QC Thresholds Explained:**
- **SNP call rate > 98% (--geno 0.02):** Ensures high-quality genotype data per marker
- **MAF > 1% (--maf 0.01):** Removes very rare variants that impute poorly and may be errors
- **HWE p > 5×10⁻⁶ (--hwe 0.000005):** Removes SNPs with significant genotyping errors (apply to controls only in case-control studies)
- **Sample call rate > 99% (--mind 0.01):** Removes poorly genotyped individuals
- **Relatedness check (--genome):** Identifies contaminated samples and cryptic relatedness
- **Sex check (--check-sex):** Identifies sample mix-ups and sex chromosome abnormalities

---

## Step 4: Strand Alignment and Allele Matching

Ensure your data matches the reference panel's strand orientation and allele coding:

```bash
# Process each chromosome
for chr in {1..22}; do
    # Normalize variants (split multi-allelic, left-align indels)
    bcftools annotate -x ID \
        -I +'%CHROM:%POS:%REF:%ALT' \
        your_data_chr${chr}_qc.vcf.gz | \
    bcftools norm -Oz -o your_data_chr${chr}_normalized.vcf.gz
    
    tabix -p vcf your_data_chr${chr}_normalized.vcf.gz
    
    # Use conform-gt to match reference panel
    java -jar conform-gt.jar \
        ref=ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        gt=your_data_chr${chr}_normalized.vcf.gz \
        match=POS \
        out=your_data_chr${chr}_conformed
done
```

**What this does:**
- Removes strand ambiguous SNPs (A/T, C/G)
- Flips alleles to match reference
- Removes SNPs not in reference panel

---

## Step 5: Pre-phasing Your Data

Phase your genotypes before imputation using either SHAPEIT4 or Eagle2:

### Option A: SHAPEIT4 (Recommended)

```bash
for chr in {1..22}; do
    shapeit4 \
        --input your_data_chr${chr}_conformed.vcf.gz \
        --map genetic_maps/chr${chr}.b37.gmap.gz \
        --region ${chr} \
        --output your_data_chr${chr}_phased.vcf.gz \
        --thread 8
done
```

### Option B: Eagle2

```bash
for chr in {1..22}; do
    eagle \
        --vcf your_data_chr${chr}_conformed.vcf.gz \
        --geneticMapFile genetic_maps/chr${chr}.b37.gmap.gz \
        --outPrefix your_data_chr${chr}_phased \
        --numThreads 8
    
    # Eagle outputs uncompressed VCF, compress it
    bgzip your_data_chr${chr}_phased.vcf
done
```

**Phasing time:** 30 minutes to 4 hours per chromosome (depends on sample size)

---

## Step 6: Imputation with Minimac4 or IMPUTE5

Perform the actual imputation using your phased genotypes and the reference panel:

### Option A: Minimac4

```bash
for chr in {1..22}; do
    minimac4 \
        --refHaps ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --haps your_data_chr${chr}_phased.vcf.gz \
        --prefix your_data_chr${chr}_imputed \
        --format GT,DS,GP \
        --allTypedSites \
        --minRatio 0.01
done
```

**Output formats:**
- **GT:** Genotypes (0/0, 0/1, 1/1)
- **DS:** Dosage (0-2 scale)
- **GP:** Genotype probabilities

### Option B: IMPUTE5 (Faster for large datasets)

```bash
for chr in {1..22}; do
    impute5 \
        --h ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --g your_data_chr${chr}_phased.vcf.gz \
        --m genetic_maps/chr${chr}.b37.gmap.gz \
        --o your_data_chr${chr}_imputed.vcf.gz \
        --r ${chr} \
        --threads 8
done
```

**Imputation time:** 1-8 hours per chromosome (depends on sample size)

---

## Step 7: Post-Imputation Quality Control

Filter imputed variants by quality metrics and perform additional sample-level QC.

### Step 7A: Variant-Level Post-Imputation QC

```bash
for chr in {1..22}; do
    # Filter by imputation quality score (R² > 0.7)
    # This is more stringent than typical thresholds (0.3-0.5)
    bcftools filter \
        -i 'R2>0.7' \
        your_data_chr${chr}_imputed.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_r2filtered.vcf.gz
    
    # Apply MAF filter (MAF > 1%)
    bcftools view \
        -i 'MAF>0.01' \
        your_data_chr${chr}_imputed_r2filtered.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_maf.vcf.gz
    
    # Convert to PLINK format for HWE testing
    plink --vcf your_data_chr${chr}_imputed_maf.vcf.gz \
          --make-bed \
          --out your_data_chr${chr}_imputed_plink
    
    # Apply HWE filter (p > 5×10⁻⁶ in controls)
    # If case-control data, use --filter-controls
    plink --bfile your_data_chr${chr}_imputed_plink \
          --hwe 0.000005 \
          --make-bed \
          --out your_data_chr${chr}_imputed_filtered
    
    # Convert back to VCF
    plink --bfile your_data_chr${chr}_imputed_filtered \
          --recode vcf \
          --out your_data_chr${chr}_imputed_final
    
    bgzip your_data_chr${chr}_imputed_final.vcf
    tabix -p vcf your_data_chr${chr}_imputed_final.vcf.gz
    
    echo "Chromosome ${chr} post-imputation QC complete"
done
```

### Step 7B: Merge Chromosomes for Sample QC

```bash
# Merge all autosomes for sample-level QC
for chr in {1..22}; do
    echo "your_data_chr${chr}_imputed_filtered"
done > merge_list.txt

plink --merge-list merge_list.txt \
      --make-bed \
      --out your_data_imputed_all_autosomes
```

### Step 7C: LD Pruning for Sample QC

Perform LD pruning to get independent SNPs for sample QC:

```bash
# Prune to independent SNPs
# Window size: 1500 kb, step: 150 SNPs, r²: 0.2
plink --bfile your_data_imputed_all_autosomes \
      --indep-pairwise 1500 150 0.2 \
      --out pruned_snps

# Extract pruned SNPs
plink --bfile your_data_imputed_all_autosomes \
      --extract pruned_snps.prune.in \
      --make-bed \
      --out your_data_imputed_pruned
```

### Step 7D: Sample-Level Post-Imputation QC

```bash
# 1. Recheck sex concordance on imputed data
plink --bfile your_data_imputed_pruned \
      --check-sex \
      --out sex_check_postimpute

# Identify sex mismatches (reported sex ≠ genetic sex)
awk 'NR>1 && $3 != 0 && $4 != $3 {print $1, $2}' sex_check_postimpute.sexcheck > sex_mismatch.txt

# 2. Check for close relatives (up to 4th degree)
# Use KING-robust for relatedness estimation
plink --bfile your_data_imputed_pruned \
      --make-king-table \
      --out relatedness_king

# Remove one from each pair with kinship > 0.088 (4th degree relatives)
plink --bfile your_data_imputed_pruned \
      --king-cutoff relatedness_king 0.088 \
      --out related_samples

# This creates related_samples.king.cutoff.out.id with samples to keep

# 3. Check heterozygosity
plink --bfile your_data_imputed_pruned \
      --het \
      --out heterozygosity

# Calculate mean and SD, remove samples beyond mean ± 3 SD
Rscript - <<'EOF'
het <- read.table("heterozygosity.het", header=TRUE)
het$HET_RATE <- (het$N.NM. - het$O.HOM.) / het$N.NM.
mean_het <- mean(het$HET_RATE)
sd_het <- sd(het$HET_RATE)
het_outliers <- het[het$HET_RATE < (mean_het - 3*sd_het) | 
                    het$HET_RATE > (mean_het + 3*sd_het), c("FID", "IID")]
write.table(het_outliers, "het_outliers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
EOF

# 4. Ancestry outlier detection (Principal Component Analysis)
# This assumes you want to keep European ancestry samples
# First, download 1000 Genomes reference populations for PCA

# Merge your data with 1000G reference (use common SNPs)
# Perform PCA
plink --bfile your_data_imputed_pruned \
      --pca 10 \
      --out pca_results

# Plot PCs and identify outliers (typically PC1 vs PC2)
# Samples clustering away from European reference = non-European ancestry
Rscript - <<'EOF'
pca <- read.table("pca_results.eigenvec", header=FALSE)
# Define European reference cluster boundaries from 1000G EUR populations
# Identify samples outside 4 SD from EUR cluster mean
# This is a simplified example - adjust based on your reference populations
# Save outliers to ancestry_outliers.txt
EOF

# 5. Combine all samples to remove
cat sex_mismatch.txt > samples_to_remove_final.txt
# Add samples NOT in related_samples.king.cutoff.out.id
cat het_outliers.txt >> samples_to_remove_final.txt
cat ancestry_outliers.txt >> samples_to_remove_final.txt

# Remove duplicates
sort samples_to_remove_final.txt | uniq > samples_to_remove_unique.txt

# Apply final sample exclusions to all chromosomes
for chr in {1..22}; do
    plink --bfile your_data_chr${chr}_imputed_filtered \
          --remove samples_to_remove_unique.txt \
          --make-bed \
          --out your_data_chr${chr}_final
    
    # Convert to VCF
    plink --bfile your_data_chr${chr}_final \
          --recode vcf \
          --out your_data_chr${chr}_final
    
    bgzip your_data_chr${chr}_final.vcf
    tabix -p vcf your_data_chr${chr}_final.vcf.gz
done
```

### Step 7E: Final QC Summary Report

```bash
# Generate comprehensive QC report
cat > final_qc_report.sh << 'SCRIPT'
#!/bin/bash

echo "=== POST-IMPUTATION QC SUMMARY ===" > post_imputation_qc_report.txt
echo "" >> post_imputation_qc_report.txt

# Count samples removed at each step
sex_mismatch=$(wc -l < sex_mismatch.txt)
related=$(wc -l < samples_to_remove_final.txt | grep -c "related")
het_outliers=$(wc -l < het_outliers.txt)
ancestry_outliers=$(wc -l < ancestry_outliers.txt)
total_removed=$(wc -l < samples_to_remove_unique.txt)

echo "Samples removed due to:" >> post_imputation_qc_report.txt
echo "  - Sex mismatch: ${sex_mismatch}" >> post_imputation_qc_report.txt
echo "  - Relatedness (>4th degree): ${related}" >> post_imputation_qc_report.txt
echo "  - Heterozygosity outliers (±3SD): ${het_outliers}" >> post_imputation_qc_report.txt
echo "  - Non-European ancestry: ${ancestry_outliers}" >> post_imputation_qc_report.txt
echo "  - Total unique samples removed: ${total_removed}" >> post_imputation_qc_report.txt
echo "" >> post_imputation_qc_report.txt

# Count final variants per chromosome
echo "Final variant counts per chromosome:" >> post_imputation_qc_report.txt
for chr in {1..22}; do
    nvars=$(wc -l < your_data_chr${chr}_final.bim)
    echo "  Chr ${chr}: ${nvars}" >> post_imputation_qc_report.txt
done

# Total autosomal variants
total_vars=$(cat your_data_chr{1..22}_final.bim | wc -l)
echo "" >> post_imputation_qc_report.txt
echo "Total autosomal variants: ${total_vars}" >> post_imputation_qc_report.txt

# Final sample count
final_samples=$(wc -l < your_data_chr1_final.fam)
echo "Final sample count: ${final_samples}" >> post_imputation_qc_report.txt

cat post_imputation_qc_report.txt
SCRIPT

chmod +x final_qc_report.sh
./final_qc_report.sh
```

**Post-Imputation QC Parameters Summary:**

**Variant-level filters:**
- **R² > 0.7:** High-quality imputation only (more stringent than typical 0.3-0.5)
- **MAF > 1%:** Consistent with pre-imputation filtering
- **HWE p > 5×10⁻⁶:** Applied in controls for case-control studies

**Sample-level filters (on LD-pruned, independent SNPs):**
- **Sex concordance:** Remove samples where genetic sex ≠ reported sex
- **Relatedness:** Remove one from pairs with kinship > 0.088 (up to 4th degree relatives)
- **Heterozygosity:** Remove samples with het rate beyond mean ± 3 SD
- **Ancestry:** Remove non-European samples (if focusing on European ancestry)

**Expected final dataset:**
- ~7-8 million high-quality autosomal variants
- Sample count depends on initial cohort and QC exclusions

---

## Step 8: Merge Chromosomes (Optional)

Combine all chromosomes into a single file:

```bash
# Create list of all chromosome files
ls your_data_chr*_imputed_final.vcf.gz > vcf_list.txt

# Concatenate all chromosomes
bcftools concat \
    -f vcf_list.txt \
    -Oz -o your_data_all_chromosomes_imputed.vcf.gz

# Index the merged file
tabix -p vcf your_data_all_chromosomes_imputed.vcf.gz
```

**Note:** Merged file can be very large (50-200 GB). Consider keeping chromosome-specific files for analysis.

---

## Step 9: Liftover from GRCh37 to GRCh38

Since 1000 Genomes Phase 3 and your input data are on GRCh37/hg19, you'll need to convert (liftover) to GRCh38/hg38.

### Required Tools and Files

**Software:**
- **Picard Tools** - CrossMap or LiftoverVcf
- **bcftools** - Post-liftover processing
- **UCSC liftOver** - Alternative tool

**Chain file:**
```bash
# Download the hg19 to hg38 chain file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
```

**GRCh38 Reference Genome:**
```bash
# Download GRCh38 reference (needed for validation)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Index the reference
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

### Method 1: Picard LiftoverVcf (Recommended)

```bash
# Liftover each chromosome
for chr in {1..22} X; do
    java -Xmx8g -jar picard.jar LiftoverVcf \
        I=your_data_chr${chr}_imputed_final.vcf.gz \
        O=your_data_chr${chr}_imputed_final_hg38.vcf.gz \
        CHAIN=hg19ToHg38.over.chain \
        REJECT=your_data_chr${chr}_rejected_variants.vcf.gz \
        R=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        WARN_ON_MISSING_CONTIG=true
    
    # Index the lifted file
    tabix -p vcf your_data_chr${chr}_imputed_final_hg38.vcf.gz
done
```

**What this does:**
- Converts coordinates from GRCh37 to GRCh38
- Validates against GRCh38 reference sequence
- Outputs rejected variants that couldn't be lifted (typically <1%)

### Method 2: CrossMap (Alternative)

```bash
# Install CrossMap
pip install CrossMap --break-system-packages

# Liftover each chromosome
for chr in {1..22} X; do
    CrossMap.py vcf \
        hg19ToHg38.over.chain \
        your_data_chr${chr}_imputed_final.vcf.gz \
        GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        your_data_chr${chr}_imputed_final_hg38.vcf
    
    # Compress and index
    bgzip your_data_chr${chr}_imputed_final_hg38.vcf
    tabix -p vcf your_data_chr${chr}_imputed_final_hg38.vcf.gz
done
```

### Post-Liftover Quality Control

After liftover, perform critical QC steps:

```bash
for chr in {1..22} X; do
    echo "Processing chromosome ${chr}..."
    
    # 1. Fix chromosome names (ensure they match "chr1" or "1" format consistently)
    bcftools annotate --rename-chrs chr_name_conv.txt \
        your_data_chr${chr}_imputed_final_hg38.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_final_hg38_fixed.vcf.gz
    
    # 2. Sort the VCF (liftover can sometimes disorder variants)
    bcftools sort \
        your_data_chr${chr}_imputed_final_hg38_fixed.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_final_hg38_sorted.vcf.gz
    
    # 3. Remove duplicate positions (can occur during liftover)
    bcftools norm -d all \
        your_data_chr${chr}_imputed_final_hg38_sorted.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_final_hg38_clean.vcf.gz
    
    # 4. Left-align and normalize indels
    bcftools norm -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        your_data_chr${chr}_imputed_final_hg38_clean.vcf.gz \
        -Oz -o your_data_chr${chr}_imputed_final_hg38_normalized.vcf.gz
    
    # Index final file
    tabix -p vcf your_data_chr${chr}_imputed_final_hg38_normalized.vcf.gz
    
    echo "Chromosome ${chr} complete."
done
```

**Create chromosome name conversion file (if needed):**
```bash
# If your files use "1" and GRCh38 expects "chr1", create chr_name_conv.txt:
for i in {1..22} X Y; do
    echo -e "${i}\tchr${i}"
done > chr_name_conv.txt
```

### Validate Liftover Results

Check that liftover was successful:

```bash
for chr in {1..22} X; do
    echo "=== Chromosome ${chr} Liftover Summary ==="
    
    # Count variants before liftover
    before=$(bcftools view -H your_data_chr${chr}_imputed_final.vcf.gz | wc -l)
    
    # Count variants after liftover
    after=$(bcftools view -H your_data_chr${chr}_imputed_final_hg38_normalized.vcf.gz | wc -l)
    
    # Count rejected variants
    rejected=$(bcftools view -H your_data_chr${chr}_rejected_variants.vcf.gz 2>/dev/null | wc -l || echo 0)
    
    # Calculate success rate
    success_rate=$(echo "scale=2; ($after / $before) * 100" | bc)
    
    echo "Before liftover: ${before}"
    echo "After liftover: ${after}"
    echo "Rejected: ${rejected}"
    echo "Success rate: ${success_rate}%"
    echo ""
done
```

**Expected results:**
- Success rate should be >99%
- Rejected variants are typically:
  - In regions that differ significantly between builds
  - In gap regions
  - Multi-allelic sites that don't map cleanly

### Common Liftover Issues and Solutions

#### Issue 1: Chromosome Naming Mismatch
**Problem:** GRCh37 uses "1", "2", etc., while GRCh38 may use "chr1", "chr2"

**Solution:**
```bash
# Add "chr" prefix
bcftools annotate --rename-chrs <(seq 1 22; echo X; echo Y | awk '{print $1"\tchr"$1}') \
    input.vcf.gz -Oz -o output.vcf.gz
```

#### Issue 2: Duplicate Positions After Liftover
**Problem:** Multiple GRCh37 positions map to same GRCh38 position

**Solution:** Use `bcftools norm -d all` (included in QC steps above)

#### Issue 3: Strand Flips
**Problem:** Reference allele doesn't match GRCh38 reference

**Solution:**
```bash
# Fix reference alleles
bcftools +fix-ref your_data.vcf.gz -- -f GRCh38_reference.fna
```

#### Issue 4: High Rejection Rate (>5%)
**Problem:** Chain file issues or severe coordinate mismatches

**Solution:**
- Verify you're using the correct chain file (hg19ToHg38, not hg18ToHg19)
- Check that input data is truly GRCh37
- Consider excluding problematic regions before liftover

### Merge GRCh38 Chromosomes

After successful liftover and QC:

```bash
# Create list of GRCh38 files
ls your_data_chr*_imputed_final_hg38_normalized.vcf.gz > vcf_list_hg38.txt

# Merge all chromosomes
bcftools concat \
    -f vcf_list_hg38.txt \
    -Oz -o your_data_all_chromosomes_imputed_hg38.vcf.gz

# Index
tabix -p vcf your_data_all_chromosomes_imputed_hg38.vcf.gz
```

### Alternative: Use GRCh38 Reference Panel from the Start

**Important consideration:** If you need GRCh38 data, you might consider:

1. **1000 Genomes GRCh38 Version:**
   - Available at: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/
   - Already in GRCh38 coordinates
   - Avoids liftover entirely
   - **Recommended if starting a new project**

2. **TOPMed Reference Panel:**
   - Native GRCh38
   - Better imputation quality
   - Free via imputation server

**To use 1000 Genomes GRCh38 version:**
```bash
# Download GRCh38 version (high coverage)
for chr in {1..22} X; do
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
done
```

However, this requires your input data to also be in GRCh38, so you'd need to liftover your genotype data BEFORE imputation instead of after.

---

## Comparison: Your Pipeline vs. Collaborator's HRC Pipeline

This section helps you understand how your pipeline matches your collaborator's specifications.

### Reference Panel Differences

| Aspect | Your Pipeline (1000G Phase 3) | Collaborator's Pipeline (HRC) |
|--------|-------------------------------|-------------------------------|
| **Reference Panel** | 1000 Genomes Phase 3 | Haplotype Reference Consortium (HRC) |
| **Sample Size** | 2,504 individuals | ~32,000 individuals |
| **Variant Count** | ~84 million | ~39 million SNPs |
| **Availability** | Publicly available | Requires data access application |
| **Build** | GRCh37 (with GRCh38 liftover) | GRCh37 |

**Impact:** HRC has more samples (better for rare variant imputation) but fewer variants overall. 1000 Genomes includes more structural variants and indels. For common variants (MAF >1%), imputation quality should be comparable.

### QC Parameters - Exact Match

Your pipeline has been configured to exactly match the collaborator's QC:

#### Pre-Imputation QC (Matched)
✓ SNP call rate: >98% (--geno 0.02)  
✓ MAF: >1% (--maf 0.01)  
✓ HWE: p > 5×10⁻⁶ in controls (--hwe 0.000005)  
✓ Sample call rate: >99% (--mind 0.01)  
✓ Relatedness check: --genome for contamination  
✓ Sex check: --check-sex for concordance  

#### Post-Imputation QC (Matched)
✓ Imputation quality: R² > 0.7  
✓ MAF: >1% (--maf 0.01)  
✓ HWE: p > 5×10⁻⁶ in controls (--hardy)  
✓ LD pruning: --indep-pairwise 1500 150 0.2  
✓ Sex concordance check (post-imputation)  
✓ Relatedness: --king-cutoff 0.088 (4th degree)  
✓ Heterozygosity: mean ± 3 SD  
✓ Ancestry filtering: PCA-based (European)  

### Expected Outcomes Comparison

Based on your matching QC parameters:

| Metric | Collaborator's Results | Your Expected Results |
|--------|------------------------|----------------------|
| **Final Sample Count** | 4,738 samples | ~Similar (depends on your starting N) |
| **Final Variant Count** | 7,381,754 autosomal variants | 7-10 million variants* |
| **Pre-imputation exclusions** | 57 samples total | Similar proportion |
| **Post-imputation exclusions** | 156 samples (sex: 17, related: 32, het: 72, ancestry: 35) | Similar proportions |

*Your variant count may be slightly higher due to 1000G having more variants initially, but after R² >0.7 filtering, the difference should be modest.

### Key Considerations for Compatibility

1. **Imputation Quality Threshold (R² >0.7):**
   - This is more stringent than typical (R² >0.3-0.5)
   - Ensures only high-confidence variants in final dataset
   - Good for meta-analysis compatibility

2. **Reference Panel Choice:**
   - For meta-analysis: Ideally use same reference, but with matched QC, results should be comparable
   - Consider running a small overlap analysis if possible to verify concordance

3. **Build Consistency:**
   - Your collaborator used GRCh37
   - This pipeline provides both GRCh37 and GRCh38 (via liftover)
   - Use GRCh37 output for direct comparison unless collaborator has also lifted to GRCh38

### Recommendations for Collaboration

1. **Use GRCh37 output initially** for direct comparison with collaborator
2. **Compare overlapping variants** (if samples overlap) to validate pipeline concordance
3. **Document any deviations** from this pipeline for reproducibility
4. **Share QC reports** (Step 7E output) with collaborator
5. **Consider harmonizing** to the same reference panel for formal meta-analysis

---

## Key Quality Metrics to Check

### Pre-Imputation Metrics

1. **Call Rate:** Should be > 95% per SNP
2. **Minor Allele Frequency (MAF):** Check distribution
3. **Hardy-Weinberg Equilibrium:** No excess of deviations
4. **Number of SNPs:** Before and after QC

### Post-Imputation Metrics

1. **Imputation Quality (R² or INFO):**
   - Mean R² should be > 0.6
   - Distribution of R² scores
   
2. **Number of Variants:**
   - Genotyped: ~500K - 5M (depending on array)
   - After imputation: ~10-40M (1000G has ~84M variants)
   
3. **MAF Distribution:**
   - Compare before and after imputation
   - Check for appropriate rare variant enrichment

4. **Imputation Rate:**
   - Percentage of reference variants successfully imputed

---

## Computational Resources

### Storage Requirements

- **1000 Genomes reference:** ~250 GB
- **Your genotype data:** 1-50 GB (depends on sample size)
- **Intermediate files:** 100-500 GB
- **Final imputed data (GRCh37):** 50-200 GB
- **GRCh38 reference genome:** ~3 GB
- **Lifted over data (GRCh38):** 50-200 GB
- **Total recommended:** 1-2 TB free space

### Memory Requirements

- **Minimum:** 16 GB RAM
- **Recommended:** 32-64 GB RAM
- **For large cohorts (>10K samples):** 128 GB RAM

### Computation Time

Per chromosome (approximate, varies by sample size):

- **QC and preparation:** 10-30 minutes
- **Phasing:** 30 minutes - 4 hours
- **Imputation:** 1-8 hours
- **Post-processing:** 10-30 minutes
- **Liftover to GRCh38:** 30 minutes - 2 hours

**Total pipeline time:** 2-7 days for all chromosomes (can parallelize)

### CPU Requirements

- **Recommended:** 8-16 cores
- Most steps are parallelizable across chromosomes
- Within-chromosome parallelization available in SHAPEIT4 and IMPUTE5

---

## Troubleshooting Common Issues

### Issue 1: Memory Errors During Imputation
**Solution:** Reduce chunk size, use IMPUTE5 instead of Minimac4, or request more memory

### Issue 2: Strand Alignment Failures
**Solution:** Manually check and flip ambiguous SNPs, ensure reference genome build matches (GRCh37/hg19 for 1000G Phase 3)

### Issue 3: Low Imputation Quality
**Solution:** 
- Check if your array has good coverage of the genome
- Verify ancestry matches reference panel populations
- Consider using TOPMed reference for better performance

### Issue 4: VCF Format Errors
**Solution:** Validate VCF files with `bcftools view -H file.vcf.gz | head` and check for malformed entries

---

## Best Practices

1. **Always perform thorough QC** before imputation - garbage in, garbage out
2. **Document your QC thresholds** for reproducibility
3. **Keep intermediate files** until you verify final results
4. **Run a test chromosome first** (e.g., chr22) to validate pipeline
5. **Check imputation quality** before proceeding with analysis
6. **Use appropriate R²/INFO thresholds** for your downstream analysis
7. **Consider population stratification** - 1000G works best for populations well-represented in the reference

---

## References

1. 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526(7571), 68-74.

2. Das S, et al. (2016). Next-generation genotype imputation service and methods. Nature Genetics, 48(10), 1284-1287.

3. Loh PR, et al. (2016). Reference-based phasing using the Haplotype Reference Consortium panel. Nature Genetics, 48(11), 1443-1448.

---

## Additional Resources

- **1000 Genomes FTP:** http://ftp.1000genomes.ebi.ac.uk/
- **Michigan Imputation Server:** https://imputationserver.sph.umich.edu/ (easier alternative)
- **TOPMed Imputation Server:** https://imputation.biodatacatalyst.nhlbi.nih.gov/
- **PLINK Documentation:** https://www.cog-genomics.org/plink/
- **bcftools Manual:** http://samtools.github.io/bcftools/bcftools.html

---

## Contact and Support

For questions about this pipeline, consult your bioinformatics core or reach out to relevant software support forums.

**Pipeline Version:** 1.0  
**Last Updated:** February 27, 2026
