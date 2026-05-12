# pipeline_stage3

Stage 3 performs post-imputation QC and converts the final study outputs into PLINK2 format.

This stage consumes stage-2 imputed VCFs from `analysis/<STUDY>/stage2/`, annotates rsIDs from dbSNP GRCh38, filters variants, merges chromosomes, performs sample QC, and writes the final study-level PLINK2 dataset to `analysis/<STUDY>/stage3/`.

## 1. Stage 3 Scope

Stage 3 inputs are:
- `analysis/<STUDY>/stage2/` (Imputed VCFs + metrics)
- `analysis/<STUDY>/stage1/` (Original FAM files)
- dbSNP GRCh38 VCF + index

Stage 3 produces a finalized, analysis-ready dataset and a consolidated **Master Report** aggregating metadata from all three stages of the pipeline.

Results are organized into a strict hierarchy under `analysis/<STUDY>/stage3/`:
- **`final/`**: The definitive, QC-passed PLINK2 files.
- **`prep_chrom/`**: Intermediate per-chromosome PLINK2 files.
- **`sample_review/`**: All sample-level QC artifacts (PCA, KING, Het, Sexcheck).
- **`report/`**: Consolidated dashboards, tables, and figures.

## 2. How To Run Stage 3

Full run (HWE and Ancestry Exclusion are ENABLED by default):
```bash
sbatch src/006_stage3.sh
```

To run a specific study:
```bash
sbatch src/006_stage3.sh --study Glbd_01
```

## 3. Ordered Execution Inside Stage 3

Stage 3 is defined in `pipeline_stage3/main.nf` and uses four active modules:

- `prepare_dbsnp.nf`
- `prepare_chrom.nf`
- `sample_qc.nf`
- `finalize_study.nf`

### 3.1 Step 1: Prepare dbSNP Per Chromosome

`PREP_DBSNP_CHROM` splits the dbSNP GRCh38 VCF into one chromosome-specific file per chromosome.

For each chromosome, it:

1. selects the correct NCBI `NC_` contig name
2. extracts that contig from dbSNP
3. renames the contig to UCSC style (`chr1`, `chr2`, ..., `chrX`)
4. writes an indexed per-chromosome dbSNP VCF

This makes downstream rsID annotation deterministic and chromosome-local.

### 3.2 Step 2: Read Stage-2 Imputed VCFs

The workflow discovers:

- `analysis/*/stage2/*_chr*_GxS.imputed.vcf.gz`

It parses:

- study ID
- chromosome
- VCF path
- VCF index path

The resulting stage-2 study/chromosome objects are then joined to the dbSNP chromosome objects created in Step 1.

### 3.3 Step 3: Annotate rsIDs

`PREP_CHROM` first overlays dbSNP onto each imputed VCF:

1. `bcftools annotate -a dbsnp_chr<CHR>.vcf.gz -c ID`
2. write an annotated chromosome VCF

At this point, variants with an exact dbSNP match receive their rsID.

### 3.4 Step 4: Normalize Variant IDs And Write A Mapping File

After dbSNP annotation, `normalize_variant_ids.py` standardizes the chromosome-level variant IDs.

The normalization policy is:

1. keep the rsID when dbSNP provides one
2. otherwise use `chr:pos:REF:ALT`

This step also writes:

- `<STUDY>_chr<CHR>.variant_id_map.tsv.gz`

That mapping file is the audit trail between the original stage-2 ID, the dbSNP annotation result, and the final stage-3 ID.

### 3.5 Step 5: Filter By Imputation Quality And MAF

Still within `PREP_CHROM`, the annotated and normalized VCF is filtered with:

- `INFO/R2 >= params.min_r2`
- `INFO/MAF >= params.maf`

By default, stage 3 uses:

- `min_r2: 0.3`
- `maf: 0.01`

The output of this step is the chromosome-level filtered VCF that feeds PLINK2 conversion.

### 3.6 Step 6: Convert Each Chromosome To PLINK2 Format

`PREP_CHROM` imports the filtered imputed VCF with `plink2` using:

- `dosage=HDS`
- `--make-pgen`
- `--sort-vars`

ChrX receives additional handling with:

- `--split-par b38`
- `--lax-chrx-import`

This step produces:

- `<STUDY>_chr<CHR>.pgen`
- `<STUDY>_chr<CHR>.pvar`
- `<STUDY>_chr<CHR>.psam`

### 3.7 Step 7: Apply HWE Filtering

HWE filtering is enabled by default in stage 3.

Current behavior:

- autosomes: HWE is applied with `plink2 --hwe`
- chrX: HWE is not applied

Default parameters:

- `run_hwe: true`
- `hwe_p: 0.000005`
- `hwe_k: 0`

The per-chromosome variant QC table records:

- input variants
- rsID-assigned variants
- fallback-ID variants
- post-R2/MAF variants
- final variants after HWE

### 3.8 Step 8: Merge All Chromosomes Per Study

`SAMPLE_QC` begins by merging the chromosome-level PLINK2 files into one study-level dataset:

- `<STUDY>_allchr.pgen`
- `<STUDY>_allchr.pvar`
- `<STUDY>_allchr.psam`

This merge is the working study-level dataset used for all downstream sample QC.

### 3.9 Step 9: Update Sex From Stage 1

Stage 3 uses the stage-1 `.fam` file as the authoritative sample-sex source.

It extracts:

- `FID`
- `IID`
- sex code

and updates the merged stage-3 dataset with `plink2 --update-sex`.

### 3.10 Step 10: Build The Sample-QC Inputs

Stage 3 then creates the study-level QC assets needed for heterozygosity, relatedness, sex check, and ancestry:

1. whole-genome PLINK bed set
2. autosome-only PLINK bed set
3. LD-pruned autosome set

These are used because some sample-QC procedures are whole-genome and others are intended to be autosome-only.

### 3.11 Step 11: Run Sample QC

`SAMPLE_QC` performs four sample-QC procedures:

1. sex check
2. relatedness estimation / KING table
3. heterozygosity
4. PCA for ancestry outlier detection

The pipeline then runs `identify_sample_outliers.py` to convert the QC outputs into study-level exclusion lists:

- sex mismatches
- heterozygosity outliers
- ancestry outliers

Related sample removal is driven from the KING cutoff output.

### 3.12 Step 12: Identify Ancestry Outliers For Every Study

Stage 3 always identifies ancestry outliers.

This is mandatory.

The key design choice is that ancestry outlier exclusion is controlled separately:

- ancestry outliers are always detected
- ancestry outliers are removed only when `--exclude-ancestry-outliers` is set

This preserves a full ancestry-QC record even when the final dataset keeps those samples.

### 3.13 Step 13: Build The Final Removal List

The removal list always includes:

- sex mismatches
- heterozygosity outliers
- related samples

It conditionally includes:

- ancestry outliers when `exclude_ancestry_outliers=true`

This list is written to:

- `<STUDY>.samples_to_remove.id`

### 3.14 Step 14: Finalize The Study-Level PLINK2 Dataset

`FINALIZE_STUDY` applies the final removal list to the merged stage-3 pfile set.

If the removal list is empty, the merged dataset is copied through unchanged.

If the removal list is non-empty, the study is rewritten after `--remove`.

The final published outputs are:

- `analysis/<STUDY>/stage3/<STUDY>.pgen`
- `analysis/<STUDY>/stage3/<STUDY>.pvar`
- `analysis/<STUDY>/stage3/<STUDY>.psam`

## 4. Master Reporting System

Stage 3 introduces the **Master Report** (`master_report.py`), which provides a unified view of the entire genetics processing history for each study.

### 4.1 Data Lineage
The report tracks sample and variant counts from:
- **Stage 1**: Raw data harmonization and liftover.
- **Stage 2**: Imputation performance and variant density.
- **Stage 3**: Final variant filtering (HWE) and sample QC.

### 4.2 Visualizations
- **Ancestry Projection**: PC1 vs PC2 plot identifying ancestry outliers.
- **Imputation Metrics**: Distribution of R2 scores across chromosomes.
- **QC Metrics**: Heterozygosity and kinship distributions.

## 5. Output Conventions

### 5.1 Final Dataset (`final/`)
- `<STUDY>.pgen`, `<STUDY>.pvar`, `<STUDY>.psam`

### 5.2 Intermediate Process Assets
- **`prep_chrom/`**: Per-chromosome PGENs after R2/MAF/HWE filtering.
- **`sample_review/`**:
    - `*.eigenvec` / `*.eigenval`: PCA results.
    - `*.king.kin0`: Kinship coefficients.
    - `*.het`: Heterozygosity statistics.
    - `*.sexcheck`: PLINK sex concordance results.

### 5.3 Reporting Assets (`report/`)
- **`report-master.html`**: The unified study dashboard.
- **`figures/`**: Ancestry PCA plots and QC distributions.
- **`tables/`**: Consolidated variant and sample metric TSVs.
- **`flags/`**: Lists of samples flagged for sex mismatches, relatedness, or ancestry.
- **`manifests/`**: Final variant ID mappings and sample removal lists.

