# pipeline_stage3

Stage 3 performs post-imputation QC and converts the final study outputs into PLINK2 format.

This stage consumes stage-2 imputed VCFs from `analysis/<STUDY>/stage2/`, annotates rsIDs from dbSNP GRCh38, filters variants, merges chromosomes, performs sample QC, and writes the final study-level PLINK2 dataset to `analysis/<STUDY>/stage3/final/`.

## 1. Stage 3 Scope

Stage 3 inputs are:
- `analysis/<STUDY>/stage2/` (Imputed VCFs + metrics)
- `analysis/<STUDY>/stage1/` (Original FAM files)
- dbSNP GRCh38 VCF + index

Stage 3 produces a finalized, analysis-ready dataset and a Stage 3-only report. Cross-stage master reports are generated separately with `src/007_report.sh`.

Results are organized into a strict hierarchy under `analysis/<STUDY>/stage3/`:
- **`final/`**: The definitive, QC-passed PLINK2 files.
- **`prep_chrom/`**: Optional intermediate per-chromosome PLINK2 files when `--publish-intermediate-plink` is enabled.
- **`sample_review/`**: All sample-level QC artifacts (PCA, KING, Het, Sexcheck).
- **`report/`**: Consolidated dashboards, tables, and figures.

## 2. How To Run Stage 3

Full run (HWE filtering and ancestry outlier exclusion are both enabled by default):
```bash
sbatch src/006_stage3.sh
```

To run a specific study:
```bash
sbatch src/006_stage3.sh --study Glbd_01
```

To retain ancestry outliers in the final dataset (detection still runs):
```bash
sbatch src/006_stage3.sh --no-exclude-ancestry-outliers
```

Large intermediate PLINK/PGEN/BED files are not copied into `analysis/` by default. To publish them for debugging or handoff:
```bash
sbatch src/006_stage3.sh --publish-intermediate-plink
```

## 3. Ordered Execution Inside Stage 3

Stage 3 is defined in `pipeline_stage3/main.nf` and uses four active modules:

- `prepare_dbsnp.nf`
- `prepare_chrom.nf`
- `sample_review.nf` (`MERGE_STUDY`, `SEX_CHECK`, `PRUNE_AUTOSOMES`, `SAMPLE_QC`, `SAMPLE_REVIEW_SUMMARY`)
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
2. stream the annotated VCF into variant-ID normalization

At this point, variants with an exact dbSNP match receive their rsID.

### 3.4 Step 4: Normalize Variant IDs And Write A Mapping File

After dbSNP annotation, `normalize_variant_ids.py` standardizes the chromosome-level variant IDs.

The normalization policy is:

1. keep the rsID when dbSNP provides one
2. otherwise use `chr:pos:REF:ALT`

This step also writes:

- `<STUDY>_chr<CHR>.variant_id_map.tsv.gz`

That mapping file is the audit trail between the original stage-2 ID, the dbSNP annotation result, and the final stage-3 ID.

Annotation, ID normalization, and R2/MAF filtering are streamed together so stage 3 does not write full annotated or normalized temporary VCFs.

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

`MERGE_STUDY` begins sample review by merging the chromosome-level PLINK2 files into one study-level dataset:

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

and writes a sex-update file. `SEX_CHECK` applies it to the chrX QC bed set for sex check, and `FINALIZE_STUDY` applies it when writing the final stage-3 PLINK2 dataset.

### 3.10 Step 10: Build The Sample-QC Inputs

The split sample-review processes then create the study-level QC assets needed for heterozygosity, relatedness, sex check, and ancestry:

1. `SEX_CHECK`: chromosome-X PLINK bed set for sex check when chrX is present
2. `PRUNE_AUTOSOMES`: LD-pruned autosome set
3. `SAMPLE_QC`: relatedness, heterozygosity, and PCA from the pruned set
4. `SAMPLE_REVIEW_SUMMARY`: final sample-exclusion manifests

The full merged PLINK2 set is not converted to unpruned PLINK bed format; sample QC only materializes the smaller datasets required by each QC procedure.

### 3.11 Step 11: Run Sample QC

The sample-review subprocesses perform four sample-QC procedures:

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

Stage 3 always identifies ancestry outliers regardless of whether they will be excluded.

The detection and exclusion steps are deliberately separated:

- ancestry outlier detection always runs and produces the PCA plot and outlier list
- ancestry outlier removal is enabled by default (`exclude_ancestry_outliers: true`); pass `--no-exclude-ancestry-outliers` to retain all samples in the final dataset

This design preserves a complete ancestry-QC record even in analyses that intentionally retain the full sample set.

### 3.13 Step 13: Build The Final Removal List

The removal list always includes:

- sex mismatches
- heterozygosity outliers
- related samples

It also includes ancestry outliers by default. Pass `--no-exclude-ancestry-outliers` to omit them:

- ancestry outliers included unless `--no-exclude-ancestry-outliers` is specified

This list is written to:

- `<STUDY>.samples_to_remove.id`

### 3.14 Step 14: Finalize The Study-Level PLINK2 Dataset

`FINALIZE_STUDY` applies the final removal list to the merged stage-3 pfile set.

`FINALIZE_STUDY` always applies the stage-1 sex update while writing the final dataset. If the removal list is non-empty, the study is also rewritten after `--remove`.

The final published outputs are:

- `analysis/<STUDY>/stage3/final/<STUDY>.pgen`
- `analysis/<STUDY>/stage3/final/<STUDY>.pvar`
- `analysis/<STUDY>/stage3/final/<STUDY>.psam`

## 4. Stage 3 Reporting

Stage 3 creates a Stage 3-only report at `analysis/<STUDY>/stage3/report/report-stage3.html`.
This report covers post-imputation filtering and sample review only.
Cross-stage master reports are generated separately with `src/007_report.sh`.

### 4.1 Report Contents
The Stage 3 report tracks:
- Variant ID annotation and fallback ID counts.
- R2/MAF variant filtering counts.
- HWE variant filtering counts.
- Sex, relatedness, heterozygosity, and ancestry sample-review counts.

### 4.2 Visualizations
- **Filtering Counts**: Variant and sample counts before and after Stage 3 filtering.
- **Ancestry Projection**: PC1 vs PC2 plot identifying ancestry outliers.

## 5. Output Conventions

### 5.1 Final Dataset (`final/`)
- `<STUDY>.pgen`, `<STUDY>.pvar`, `<STUDY>.psam`

### 5.2 Intermediate Process Assets
- **`prep_chrom/`**: Per-chromosome PGENs after R2/MAF/HWE filtering when `--publish-intermediate-plink` is enabled. By default, these files remain in the Nextflow work cache.
- **`sample_review/`**:
    - `*.eigenvec` / `*.eigenval`: PCA results.
    - `*.kin0`: Kinship coefficients.
    - `*.het`: Heterozygosity statistics.
    - `*.sexcheck`: PLINK sex concordance results.
    - Intermediate merged, chrX BED, and pruned BED files only when `--publish-intermediate-plink` is enabled.

### 5.3 Reporting Assets (`report/`)
- **`report-stage3.html`**: The Stage 3-only study dashboard.
- **`figures/`**: Ancestry PCA plots and QC distributions.
- **`tables/`**: Consolidated variant and sample metric TSVs.
- **`flags/`**: Lists of samples flagged for sex mismatches, relatedness, or ancestry.
- **`manifests/`**: Final variant ID mappings and sample removal lists.
