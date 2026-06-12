# Pipeline Overview

The EPIC genetics pipeline runs in four sequential stages, each submitted as a Slurm batch job. All stages read configuration from the `.env` file and write outputs to `${SCRATCH}/${SCRATCH_DATE}/`.

## Pipeline Stages

| Stage | Purpose | Submission |
| --- | --- | --- |
| [Stage 1](stage1.md) | Study harmonisation and liftover to GRCh38 | `sbatch src/004_stage1.sh` |
| [Stage 2](stage2.md) | Phasing (SHAPEIT5) and imputation (Minimac4) against 1000G GRCh38 | `sbatch src/005_stage2.sh` |
| [Stage 3](stage3.md) | Post-imputation QC, rsID annotation, PLINK2 conversion, sample QC | `sbatch src/006_stage3.sh` |
| [Stage 4](stage4.md) | Cross-stage master reports and final study archive | `bash src/007_stage4.sh` |

## Key Filters

| Filter | Stage | Threshold |
| --- | --- | --- |
| Variant MAF (pre-phasing) | Stage 2 | MAF ≥ 0.005 |
| Imputation quality | Stage 2/3 | R² ≥ 0.3 |
| Minor allele frequency | Stage 3 | MAF ≥ 0.01 |
| Hardy-Weinberg equilibrium | Stage 3 | p ≥ 0.000005 (autosomes only) |
| Heterozygosity outliers | Stage 3 | > 3.0 SD |
| Relatedness (KING) | Stage 3 | kinship ≥ 0.0884 (reported; removed only with `--exclude-related`) |
| Ancestry outliers | Stage 3 | z-score ≥ 6.0 on 10 PCs (reported; removed only with `--exclude-ancestry-outliers`) |

## Running the Full Pipeline

### 1. Configure the environment

Edit `.env` and set `SCRATCH_DATE` to your run label (e.g. `2025-06`). All outputs will be written to `${SCRATCH}/${SCRATCH_DATE}/`.

```bash
bash src/000_env.sh      # validate .env variables
sbatch src/000_tools.sh  # install pipeline dependencies (bcftools, PLINK2, etc.)
```

### 2. Prepare input data

```bash
sbatch src/001_data-genetics.sh   # copy raw genotype arrays from archive to scratch
bash src/002_data-reference.sh    # download 1000G GRCh38 reference for imputation
Rscript src/003_data-epic.R       # prepare EPIC phenotype and sample manifest files
```

### 3. Run the pipeline stages

Each stage must complete before the next is submitted. Monitor with `squeue -u $USER`.

```bash
sbatch src/004_stage1.sh   # harmonisation + liftover
sbatch src/005_stage2.sh   # phasing + imputation
sbatch src/006_stage3.sh   # post-imputation QC + PLINK2 conversion
bash src/007_stage4.sh     # reports + final archive
```

### 4. Update documentation

After all stages complete, regenerate the sample overlap plot and update summary tables:

```bash
bash src/008_documentation.sh
```

Then commit the updated `docs/` and `README.md`:

```bash
git add docs/ README.md
git commit -m "Update pipeline outputs and documentation"
git push
```

---

## Testing with a Single Study (Glbd_01)

Glbd_01 is the smallest study (N = 114) and is the recommended study for testing pipeline changes end-to-end.

### Configure `.env` for a single-study test run

Set a fresh `SCRATCH_DATE` so the test run is isolated from production outputs:

```bash
SCRATCH_DATE="test-Glbd_01"
```

### Run stages with `--study` filter

Each submission script accepts a `--study` argument that restricts the Nextflow pipeline to a single study:

```bash
sbatch src/004_stage1.sh --study Glbd_01
sbatch src/005_stage2.sh --study Glbd_01
sbatch src/006_stage3.sh --study Glbd_01
bash src/007_stage4.sh   --study Glbd_01
```

### Expected outputs

After a successful test run, check:

```
${SCRATCH}/${SCRATCH_DATE}/studies/Glbd_01/stage3/
  ├── Glbd_01.pgen / .psam / .pvar.zst   # PLINK2 final dataset
  ├── sample_qc/                          # heterozygosity, relatedness, ancestry outlier reports
  └── variants_filtered.log              # variant counts after each filter

${SCRATCH}/${SCRATCH_DATE}/final/
  ├── report/Glbd_01.master-report.html
  └── summaries/Glbd_01.summary.md
```

Run `008_documentation.sh` as normal — it will copy only the files that exist, so a single-study test will produce a Glbd_01-only report without affecting other study outputs.
