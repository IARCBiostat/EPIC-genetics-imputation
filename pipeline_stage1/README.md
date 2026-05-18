# pipeline_stage1

Stage 1 converts each raw EPIC genetics dataset into a harmonized hg38 PLINK handoff ready for stage 2.

The stage-1 workflow standardizes all studies into a common coordinate system and genome build.

## 1. Stage 1 Scope

Stage 1 starts from study-specific raw PLINK files and ends with:

- `analysis/<STUDY>/stage1/<STUDY>.bed`
- `analysis/<STUDY>/stage1/<STUDY>.bim`
- `analysis/<STUDY>/stage1/<STUDY>.fam`
- `analysis/<STUDY>/stage1/summary.txt`

Stage-1 temporary files are kept under:

- `pipeline_stage1/work/<STUDY>/Trace/`
- `pipeline_stage1/work/<STUDY>/Trace/QC_init/`
- `pipeline_stage1/work/<STUDY>/Trace/QC_completion/`
- `pipeline_stage1/work/<STUDY>/Trace/QC_exclusion/`
- `pipeline_stage1/work/<STUDY>/Trace/LiftOver/`

Stage 1 is responsible for:

1. copying the raw study data into a controlled work area
2. comparing the dataset to the study chip manifest
3. applying study-specific SNP and sample harmonization
4. standardizing sample IDs to EPIC IDs
5. executing the two-step SNP exclusion and completion logic carried over from the archived pipeline design
6. lifting the final post-QC PLINK dataset to GRCh38 / hg38

## 2. Required Input Layout

`--data-root` must point to the archive-style raw data root that contains both the study folders and the shared EPIC ID reference:

- `Breast/`
- `Colonrectum/`
- `Epic_Cvd/`
- `Gallbladder/`
- `Interact/`
- `Kidney/`
- `Lung/`
- `Lymphoma/`
- `Ovary/`
- `Pancreas/`
- `Prostate/`
- `Stomach/`
- `Uadt/`
- `Reference/Epic/Subj_Id_2015.txt`

The stage-1 scripts do not operate directly on the flattened sync layout in `data/genetics/`.

## 3. How To Run Stage 1

The normal batch entrypoint is:

```bash
sbatch src/004_stage1.sh
```

To run a single study:

```bash
STAGE1_SCRIPTS=process_brea_01_erneg.py sbatch src/004_stage1.sh
```

To run one study script directly:

```bash
python3 pipeline_stage1/scripts/process_brea_01_erneg.py \
  --data-root /path/to/archive_style_root \
  --work-root pipeline_stage1/work \
  --plink /path/to/plink \
  --python2 /path/to/python2.7
```

Required runtime tools:

- `python3`
- `python2.7`
- `plink`
- `perl`
- Linux `x86_64` support for the bundled UCSC `liftOver` binary used by the vendored `triple-liftOver` workflow

## 4. Ordered Execution Inside Each Stage-1 Script

All bespoke stage-1 scripts follow the same numbered structure. The study-specific behavior is controlled by the configuration block at the top of each script:

- `BUILD`
- `ID_LINK_REL`
- `PRE`
- `PART1`
- `COMPLETION`
- `PART2_BUILD35_REL`

### 4.1 Step 1: Create The Study Work Area

Each script creates:

- `Trace/`
- `Trace/QC_init/`
- `Trace/QC_completion/`
- `Trace/QC_exclusion/`
- `Trace/LiftOver/`
- `analysis/<STUDY>/stage1/`

This isolates the full preprocessing history for each study and ensures that the final stage-1 handoff lands in a stable downstream analysis location.

### 4.2 Step 2: Copy The Raw PLINK Prefix Into `Trace/`

The script copies the study’s raw `.bed/.bim/.fam` prefix from the archive-style data root into `Trace/<STUDY>`.

This prevents the raw source dataset from being modified directly and guarantees that every downstream operation is acting on a local working copy.

### 4.3 Step 3: Run Initial Manifest QC

Each script runs `pipeline_stage1/bin/preProcessing.py` on the copied PLINK dataset using the study-specific chip manifest.

The initial QC produces the reference comparison files that drive the downstream SNP-exclusion logic, including:

- manifest lookup outputs
- SNPs not found in the manifest
- strand disagreement files
- chromosome/position disagreement files
- allele disagreement files
- duplicate and mitochondrial diagnostics

These files are moved into `Trace/QC_init/`.

### 4.4 Step 4: Run Optional Pre-ID Corrections

Some studies need special preprocessing before EPIC ID harmonization. Depending on the study, the script may:

- resort the input PLINK dataset
- perform SNP name linkage between BIM and manifest
- reset SNP positions
- exclude X/Y/XY/MT before sample harmonization
- convert AB allele nomenclature to explicit alleles

The study-specific `PRE` configuration field controls which of these steps are applied. The logic is self-documented in the configuration block at the top of each study script under `pipeline_stage1/scripts/`.

### 4.5 Step 5: Harmonize Sample IDs To EPIC IDs

Each script then aligns sample IDs to the EPIC master ID system using:

- the current study `.fam`
- `Reference/Epic/Subj_Id_2015.txt`
- optionally a study-specific link file

This step produces:

- `*_goodID.txt`
- `*_removeID.txt`

The script then applies the linked sample state with PLINK in this order:

1. `--update-sex`
2. `--remove`
3. `--update-ids`
4. `--make-pheno` using an empty phenotype file
5. `--indiv-sort n`

The result is the study-level EPIC-harmonized dataset:

- `Trace/<STUDY>_EPIC_all`

### 4.6 Step 6: Run SNP Exclusion Part 1

Part 1 exclusion always begins by removing SNPs not found in the manifest comparison. The scripts then diverge by study for:

- negative strand handling
- chromosome/position mismatch handling
- allele mismatch handling
- optional unknown-strand exclusion

This is the first major study-specific branch point in the pipeline. The study-specific `PART1` configuration field in each script defines which exclusion classes apply.

### 4.7 Step 7: Run The Completion Stage

The completion stage repairs metadata after the first exclusion step. Depending on the study, the script may:

- search and update chromosome codes
- search and update positions
- flip strand against the appropriate manifest orientation
- search and update alleles
- use an auxiliary loci-to-rsID reference file
- apply a build 36.2 subversion correction

This stage converts the remaining variants into the best possible harmonized pre-final state before the second exclusion pass.

### 4.8 Step 8: Rerun QC After Completion

`preProcessing.py` is rerun on the completed dataset and the outputs are stored in `Trace/QC_completion/`.

These completion-QC files are the inputs to SNP Exclusion Part 2.

### 4.9 Step 9: Run SNP Exclusion Part 2

All studies apply the same core post-completion cleanup chain:

1. misplaced mitochondrial SNP exclusion
2. duplicate type 4 exclusion
3. duplicate type 2 exclusion
4. duplicate type 3 exclusion
5. duplicate type 1 exclusion

Two studies also apply an additional build-35 exclusion file after the duplicate cleanup.

The final product of this section is the study’s post-QC stage-1 prefix:

- `Trace/<STUDY>_EXC2_all`

### 4.10 Step 10: Run Final QC On The Post-Exclusion Dataset

Each script runs `preProcessing.py` one final time on `*_EXC2_all`.

These outputs are stored in `Trace/QC_exclusion/` and serve as the final pre-liftover QC audit trail.

### 4.11 Step 11: Liftover The Final Dataset To hg38

The post-QC dataset is lifted to GRCh38 / hg38 using the vendored Perl `triple-liftOver` workflow.

The stage-1 liftover step:

1. selects the correct chain based on the study build
2. updates chromosome codes and positions
3. flips inverted SNVs identified by `triple-liftOver`
4. runs `--split-x b38 no-fail`
5. writes the final lifted handoff to `analysis/<STUDY>/stage1/`

Build handling:

- build 36 studies use `hg18ToHg38`
- build 37 studies use `hg19ToHg38`

### 4.12 Step 12: Retain Trace Artefacts

Intermediate exclusion lists, QC outputs, linkage files, and liftover artefacts are intentionally retained.

This is deliberate: stage 1 is designed to be auditable study by study.

## 5. Reporting

After all study scripts complete, `src/004_stage1.sh` runs two reporting steps:

1. **Study-level reports** — `pipeline_stage1/scripts/run_stage1_reports.py` generates per-study figures, tables, and HTML report assets under `analysis/<STUDY>/stage1/report/`. These are consumed by the cross-stage master report generator (`src/007_report.sh`).

2. **Cross-study summary** — `pipeline_stage1/scripts/summary.py` writes a markdown summary of sample and variant counts across all studies to `analysis/stage1-summary.md`.

## 6. Study-Specific Configuration Fields

The stage-1 scripts encode all study-specific logic in the configuration block at the top of each file:

| Field | Meaning |
| --- | --- |
| `BUILD` | Starting genome build of the raw study data (`36`, `37`, or `38`) |
| `ID_LINK_REL` | Study-specific sample ID linkage file, if required |
| `PRE` | Pre-ID special handling such as sorting, position reset, AB translation, or chromosome exclusion |
| `PART1` | Study-specific Part 1 SNP exclusion behaviour |
| `COMPLETION` | Study-specific metadata completion behaviour |
| `PART2_BUILD35_REL` | Optional extra build-35 exclusion list applied in Part 2 |

## 7. Studies

The following studies are processed by Stage 1:

| Script | Study ID | Description |
| --- | --- | --- |
| `process_brea_01_erneg.py` | Brea_01_Erneg | Breast cancer — ER-negative |
| `process_brea_02.py` | Brea_02 | Breast cancer |
| `process_clrt_01.py` | Clrt_01 | Colorectal cancer |
| `process_ecvd_01.py` | Ecvd_01 | Cardiovascular disease 1 |
| `process_ecvd_02.py` | Ecvd_02 | Cardiovascular disease 2 |
| `process_ecvd_03.py` | Ecvd_03 | Cardiovascular disease 3 |
| `process_glbd_01.py` | Glbd_01 | Gallbladder |
| `process_inte_01.py` | Inte_01 | InterAct 1 |
| `process_inte_02.py` | Inte_02 | InterAct 2 |
| `process_inte_03.py` | Inte_03 | InterAct 3 |
| `process_kidn_01.py` | Kidn_01 | Kidney cancer 1 |
| `process_kidn_02.py` | Kidn_02 | Kidney cancer 2 |
| `process_lung_01.py` | Lung_01 | Lung cancer |
| `process_lymp_01.py` | Lymp_01 | Lymphoma |
| `process_neuro_01.py` | Neuro_01 | Neurological |
| `process_ovar_01.py` | Ovar_01 | Ovarian cancer |
| `process_panc_01.py` | Panc_01 | Pancreatic cancer 1 |
| `process_panc_02.py` | Panc_02 | Pancreatic cancer 2 |
| `process_pros_01.py` | Pros_01 | Prostate cancer 1 |

| `process_pros_03.py` | Pros_03 | Prostate cancer 3 |
| `process_pros_04.py` | Pros_04 | Prostate cancer 4 |
| `process_stom_01.py` | Stom_01 | Stomach cancer |
| `process_uadt_01.py` | Uadt_01 | Upper aerodigestive tract cancer |

