# pipeline_stage4

Stage 4 assembles the outputs from all preceding pipeline stages into a per-study deliverable. It generates a cross-stage master HTML report integrating stage 2 and stage 3 QC summaries, then packages the finalized PLINK2 files, QC exclude lists, review files, and all HTML reports into a structured archive ready for downstream analysis and dissemination.

This stage consumes stage 3 outputs from `analysis/<STUDY>/stage3/` and stage 2 reports from `analysis/<STUDY>/stage2/report/`, and writes the final deliverable to `final/<STUDY>/`.

## 1. Stage 4 Scope

Stage 4 inputs are:

- `analysis/<STUDY>/stage3/final/` — per-chromosome PLINK2 pgen/pvar files and a single psam from stage 3 finalization
- `analysis/<STUDY>/stage3/report/flags/` — QC exclude lists (related, ancestry, HWE) written by stage 3
- `analysis/<STUDY>/stage3/sample_review/` — KING kinship table, heterozygosity, and PCA files from stage 3 sample review
- `analysis/<STUDY>/stage3/report/report-stage3.html` — the per-study stage 3 QC report
- `analysis/<STUDY>/stage2/report/report-stage2.html` — the per-study stage 2 imputation report

Stage 4 outputs are organized under `final/<STUDY>/`:

- **`<STUDY>.tar.gz`**: The deliverable archive containing all per-chromosome PLINK2 files for downstream analysis.
- **`report-master.html`**: Cross-stage master HTML report integrating stage 2 and stage 3 outputs.
- **`report-stage2.html`**: Copy of the stage 2 imputation report.
- **`report-stage3.html`**: Copy of the stage 3 QC report.
- **`review/`**: QC review files (kinship, heterozygosity, PCA, and exclude lists) for analyst inspection.

The cross-stage master reports are also published to a shared `report/` directory:

- **`report/<STUDY>.master-report.html`**: Indexed master report for cross-study comparison.

## 2. How To Run Stage 4

Full run (all studies):

```bash
bash src/007_stage4.sh
```

To run a specific study:

```bash
bash src/007_stage4.sh --study Glbd_01
```

To specify a custom destination root:

```bash
bash src/007_stage4.sh --dest-root /path/to/final
```

Stage 4 must be run after stage 3 has completed for the target study or studies. It reads stage 3 outputs directly from disk rather than from a Nextflow channel, so it does not need to run in the same Nextflow session as stage 3.

### 2.1 Key Parameters

All parameters have defaults derived from `.env`. They can also be passed as flags to `007_stage4.sh`:

| Flag | Default | Description |
| --- | --- | --- |
| `--study` | `all` | Comma-separated study IDs, or `all` |
| `--analysis-root` | `${SCRATCH_RUN}/studies` | Root containing per-study stage outputs |
| `--stage2-root` | analysis-root | Root containing stage 2 outputs (if different) |
| `--stage3-root` | analysis-root | Root containing stage 3 outputs (if different) |
| `--dest-root` | `${SCRATCH_RUN}/final` | Destination for deliverable archives |
| `--report-dir` | `dest-root/report` | Directory for master HTML copies |
| `--partition` | `low_p` | Slurm partition |
| `--no-resume` | (resume on by default) | Disable Nextflow `-resume` |

### 2.2 Resource Configuration

Stage 4 is lightweight. Both processes run with minimal resources as the bottleneck is disk I/O, not computation.

| Process | CPUs | Memory | Time |
| --- | --- | --- | --- |
| `MASTER_REPORT` | 1 | 8 GB | 4h |
| `FINALISE_STUDY` | 8 | 16 GB | 4h |

`FINALISE_STUDY` requests 8 CPUs because it uses `pigz` (parallel gzip) when available, which parallelizes compression of the deliverable tarball.

## 3. Ordered Execution Inside Stage 4

Stage 4 is defined in `pipeline_stage4/main.nf` and uses two modules:

- `modules/report.nf` — `MASTER_REPORT`
- `modules/finalise.nf` — `FINALISE_STUDY`

### 3.1 Step 1: Discover Stage 3 Outputs

The workflow uses `Channel.fromPath` glob patterns to discover all stage 3 outputs on disk, deriving the study name from the path structure:

```
${stage3_root}/<STUDY>/stage3/final/*.pgen
${stage3_root}/<STUDY>/stage3/final/*.pvar
${stage3_root}/<STUDY>/stage3/final/*_chr1.psam
${stage3_root}/<STUDY>/stage3/report/flags/*.related.exclude
${stage3_root}/<STUDY>/stage3/report/flags/*.ancestry.exclude
${stage3_root}/<STUDY>/stage3/report/flags/*.hwe.exclude
${stage3_root}/<STUDY>/stage3/sample_review/*_het.het
${stage3_root}/<STUDY>/stage3/sample_review/*_king.kin0
${stage3_root}/<STUDY>/stage3/sample_review/*_pca.eigenvec
${stage3_root}/<STUDY>/stage3/sample_review/*_pca.eigenval
${stage3_root}/<STUDY>/stage3/report/report-stage3.html
```

Only the chr1 psam is staged per study; all chromosomes carry an identical sample list so chr1 is used as the canonical representative.

### 3.2 Step 2: Generate The Master Report

`MASTER_REPORT` runs `src/misc/master_report.py` for each study.

The script reads the full report tree from `analysis_root` — it does not receive report files as Nextflow inputs — and generates a single self-contained HTML document that integrates:

- stage 2 imputation metrics (variant counts, R², empirical R², allele frequency correlation)
- stage 3 QC metrics (variant filtering, sample filtering, ancestry PCA, heterozygosity)
- cross-stage filtering step summary

The stage 3 HTML path is used as the Nextflow input only as a **cache sentinel**: Nextflow uses it to determine whether the master report is up to date. The actual report content is assembled by reading the full `analysis_root` directory tree at runtime.

The master report is published to:
- `${params.report_dir}/<STUDY>.master-report.html`

### 3.3 Step 3: Finalise The Study Deliverable

`FINALISE_STUDY` receives all stage 3 outputs for a study as a single joined tuple and assembles the deliverable. It runs with `stageInMode = 'copy'` so all input files are local copies in the Nextflow work directory rather than NFS symlinks, which ensures reliable archiving.

The process performs the following steps in order:

#### 3.3.1 Build the archive directory

All per-chromosome PGEN and PVAR files are moved (not copied) into a `<STUDY>/` subdirectory within the work directory. The single canonical PSAM is copied in as `<STUDY>/<STUDY>.psam`.

The archive contains:
```
<STUDY>/
├── <STUDY>_chr1.pgen
├── <STUDY>_chr1.pvar
├── <STUDY>_chr2.pgen
├── <STUDY>_chr2.pvar
...
├── <STUDY>_chr22.pgen
├── <STUDY>_chr22.pvar
├── <STUDY>_chrX.pgen
├── <STUDY>_chrX.pvar
└── <STUDY>.psam
```

#### 3.3.2 Compress the archive

The directory is compressed to `<STUDY>.tar.gz`. If `pigz` (parallel gzip) is available on the node, it is used with all available CPUs for faster compression. Otherwise the process falls back to standard `tar -czf`.

#### 3.3.3 Copy review files

The stage 3 QC files are copied into a `review/` subdirectory with normalized names:

| Destination | Source |
| --- | --- |
| `review/related.exclude` | `*.related.exclude` from stage 3 flags |
| `review/ancestry.exclude` | `*.ancestry.exclude` from stage 3 flags |
| `review/hwe.exclude` | `*.hwe.exclude` from stage 3 flags |
| `review/het.het` | `*_het.het` from stage 3 sample review |
| `review/king.kin0` | `*_king.kin0` from stage 3 sample review |
| `review/pca.eigenvec` | `*_pca.eigenvec` from stage 3 sample review |
| `review/pca.eigenval` | `*_pca.eigenval` from stage 3 sample review |

#### 3.3.4 Copy HTML reports

- The stage 3 HTML report is copied to `report-stage3.html`.
- The stage 2 HTML report is read from `${stage2_root}/<STUDY>/stage2/report/report-stage2.html` and copied to `report-stage2.html`. If no stage 2 report is found, a minimal placeholder HTML is written instead.
- The master report generated in step 3.2 is copied to `report-master.html`.

All outputs are published to `${params.dest_root}/<STUDY>/`.

## 4. Post-Pipeline Steps

After the Nextflow pipeline completes, `007_stage4.sh` copies the per-stage summary markdown files into `${dest_root}/summaries/`:

- `summaries/stage1-summary.md`
- `summaries/stage2-summary.md`
- `summaries/stage3-summary.md`

These summaries are read from `analysis_root` where each stage script writes them on completion.

## 5. Output Conventions

### 5.1 Deliverable Archive (`<STUDY>.tar.gz`)

The tarball contains all per-chromosome PLINK2 files required for downstream analysis. Extracting it produces a `<STUDY>/` directory with:

- 23 pairs of `.pgen` / `.pvar` files (chromosomes 1–22 and X)
- one `.psam` file (sample list, identical across all chromosomes)

The PSAM contains updated sex codes from the stage 1 `.fam` and phenotype information. The PGEN files contain dosage data at all variants that passed stage 3 R²/MAF filtering, with the sample exclusions applied by stage 3 `FINALIZE_CHROM`.

### 5.2 Review Files (`review/`)

The `review/` directory contains the stage 3 QC artifacts required for post-hoc analyst review:

| File | Contents |
| --- | --- |
| `related.exclude` | Sample IDs identified as related (KING kinship ≥ 0.0884) |
| `ancestry.exclude` | Sample IDs identified as ancestry outliers (PCA z-score ≥ 6.0 on 10 PCs) |
| `hwe.exclude` | Variant IDs failing HWE (p < 0.000005) on autosomes; for downstream use if needed |
| `het.het` | Per-sample autosomal heterozygosity statistics from the LD-pruned dataset |
| `king.kin0` | Pairwise kinship coefficients (KING format) |
| `pca.eigenvec` | Per-sample PCA scores (10 PCs) |
| `pca.eigenval` | PCA eigenvalues |

### 5.3 HTML Reports

| File | Contents |
| --- | --- |
| `report-stage2.html` | Stage 2 imputation report (phasing quality, R², empirical R², AF correlation) |
| `report-stage3.html` | Stage 3 QC report (variant filtering, sample QC, ancestry PCA, HWE counts) |
| `report-master.html` | Cross-stage master report integrating stage 2 and stage 3 summaries |

### 5.4 Configuration and Environment

Stage 4 uses a minimal conda environment (`envs/report.yml`) containing only Python ≥ 3.9 for the master report script. `FINALISE_STUDY` uses no conda environment — it relies only on standard system tools (`bash`, `tar`, `cp`, `pigz` if available).
