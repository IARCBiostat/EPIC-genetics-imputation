# pipeline_stage1

Stage 1 converts each raw EPIC genetics dataset into a harmonized hg38 PLINK handoff ready for stage 2.

This stage is implemented as one bespoke script per study in `pipeline_stage1/scripts/`. The archived 2022 scripts in `temp/archive/pipeline-2022/` are used only as templates and references; they are not called directly by the active stage-1 workflow.

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
sbatch src/003_stage1.sh
```

To run a single study:

```bash
STAGE1_SCRIPTS=process_brea_01_erneg.py sbatch src/003_stage1.sh
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

The detailed study-specific matrix appears in Section 6.

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

This is the first major study-specific branch point in the pipeline and is documented in Sections 6 and 7.

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

## 5. Study-Specific Configuration Fields

The stage-1 scripts encode all study-specific logic in the configuration block at the top of each file:

| Field | Meaning |
| --- | --- |
| `BUILD` | Starting genome build of the raw study data |
| `ID_LINK_REL` | Study-specific sample ID linkage file, if required |
| `PRE` | Pre-ID special handling such as sorting, position reset, AB translation, or chromosome exclusion |
| `PART1` | Study-specific Part 1 SNP exclusion behavior |
| `COMPLETION` | Study-specific metadata completion behavior |
| `PART2_BUILD35_REL` | Optional extra build-35 exclusion list applied in Part 2 |

## 6. Study-Specific Matrix

This table summarizes the study-specific configuration actually encoded in the bespoke scripts.

| Study | Build | ID linkage | Pre-ID step(s) | Part 1 | Completion | Part 2 special |
| --- | --- | --- | --- | --- | --- | --- |
| Brea_01_Erneg | 37 | Study link file | None | `neg=exclude_FR`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_fr_pm`, `alleles=search` | None |
| Brea_02 | 37 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_tb_pm`, `alleles=search` | None |
| Clrt_01 | 37 | Study link file | `name_linkage` | `neg=exclude_PM`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=move`, `alleles=move` | None |
| Ecvd_01 | 37 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_tb_pm`, `alleles=search` | None |
| Ecvd_02 | 37 | Study link file | `reset_positions` | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude_same`, `unknown=PM` | `chr=search`, `pos=search`, `strand=flip_tb_pm`, `alleles=search` | None |
| Ecvd_03 | 37 | Study link file | `ab_translate` | `neg=move`, `chrpos=exclude`, `alleles=move_same`, `unknown=PM` | `chr=move`, `pos=move`, `strand=flip_pm`, `alleles=move` | None |
| Glbd_01 | 37 | Study link file | None | `neg=exclude_FR`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_fr_pm`, `alleles=search` | None |
| Inte_01 | 36 | EPIC master only | None | `neg=exclude_PM`, `chrpos=move`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=move`, `alleles=move` | None |
| Inte_02 | 37 | EPIC master only | None | `neg=exclude_PM`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=move`, `alleles=move` | None |
| Inte_03 | 37 | EPIC master only | None | `neg=exclude_PM`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=move`, `alleles=search` | None |
| Kidn_01 | 36 | Study link file | None | `neg=exclude_PM`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=move`, `strand=move`, `alleles=search` | None |
| Kidn_02 | 37 | EPIC master only | None | `neg=exclude_PM`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=search`, `strand=move`, `alleles=search`, `txt_rel` | None |
| Lung_01 | 37 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_tb_pm`, `alleles=search` | None |
| Lymp_01 | 36 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_tb_pm`, `alleles=search`, `subversion_rel` | build-35 exclusion |
| Ovar_01 | 37 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_tb_pm`, `alleles=search` | None |
| Panc_01 | 36 | Study link file | None | `neg=skip`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_pm`, `alleles=search` | None |
| Panc_02 | 36 | Study link file | None | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_tb_pm`, `alleles=search`, `subversion_rel` | build-35 exclusion |
| Pros_01 | 36 | Study link file | `sort_input` | `neg=exclude_FR`, `chrpos=skip`, `alleles=exclude` | `chr=move`, `pos=move`, `strand=flip_fr_pm`, `alleles=move` | None |
| Pros_02 | 37 | Study link file | None | `neg=skip`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_pm`, `alleles=search` | None |
| Pros_03 | 37 | Study link file | `reset_positions` | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_tb_pm`, `alleles=search` | None |
| Pros_04 | 37 | Study link file | `exclude_xymt` | `neg=exclude_TB`, `chrpos=exclude`, `alleles=exclude` | `chr=search`, `pos=search`, `strand=flip_tb_pm`, `alleles=search` | None |
| Stom_01 | 36 | Study link file | None | `neg=exclude_TB`, `chrpos=move`, `alleles=exclude` | `chr=move`, `pos=search`, `strand=flip_tb_pm`, `alleles=move` | None |
| Uadt_01 | 37 | EPIC master only | `ab_translate`, `ab_txt_rel` | `neg=move`, `chrpos=exclude`, `alleles=move_all` | `chr=move`, `pos=move`, `strand=flip_pm`, `alleles=search` | None |

## 7. Current Stage-1 Summary Snapshot

The table below mirrors the current `analysis/stage1-summary.md` summary snapshot.

### 7.1 Complete Studies

| Study | Samples | Variants | Sex 1 | Sex 2 | Sex Other | Pheno 1 | Pheno 2 | Pheno Other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Brea_01_Erneg | 1,011 | 558,570 | 0 | 1,011 | 0 | 0 | 0 | 1,011 |
| Brea_02 | 7,491 | 478,002 | 0 | 7,491 | 0 | 0 | 0 | 7,491 |
| Ecvd_01 | 9,426 | 531,696 | 5,386 | 4,040 | 0 | 0 | 0 | 9,426 |
| Ecvd_02 | 8,920 | 208,221 | 5,179 | 3,741 | 0 | 0 | 0 | 8,920 |
| Ecvd_03 | 8,587 | 395,482 | 5,037 | 3,550 | 0 | 0 | 0 | 8,587 |
| Glbd_01 | 119 | 689,010 | 25 | 94 | 0 | 0 | 0 | 119 |
| Inte_01 | 9,290 | 566,795 | 3,902 | 5,388 | 0 | 0 | 0 | 9,290 |
| Inte_02 | 7,397 | 516,117 | 2,946 | 4,451 | 0 | 0 | 0 | 7,397 |
| Inte_03 | 6,328 | 509,114 | 3,153 | 3,175 | 0 | 0 | 0 | 6,328 |
| Kidn_01 | 356 | 580,897 | 195 | 161 | 0 | 0 | 0 | 356 |
| Kidn_02 | 265 | 4,146,971 | 144 | 121 | 0 | 0 | 0 | 265 |
| Lung_01 | 2,484 | 492,643 | 1,549 | 935 | 0 | 0 | 0 | 2,484 |
| Lymp_01 | 480 | 732,277 | 224 | 256 | 0 | 0 | 0 | 480 |
| Ovar_01 | 1,310 | 470,489 | 0 | 1,310 | 0 | 0 | 0 | 1,310 |
| Panc_01 | 751 | 559,366 | 380 | 371 | 0 | 0 | 0 | 751 |
| Panc_02 | 183 | 732,277 | 71 | 112 | 0 | 0 | 0 | 183 |
| Pros_01 | 856 | 567,416 | 856 | 0 | 0 | 0 | 0 | 856 |
| Pros_02 | 1,801 | 201,475 | 1,801 | 0 | 0 | 0 | 0 | 1,801 |
| Pros_03 | 1,137 | 496,063 | 1,137 | 0 | 0 | 0 | 0 | 1,137 |
| Pros_04 | 1,488 | 669,060 | 1,488 | 0 | 0 | 0 | 0 | 1,488 |
| Stom_01 | 317 | 654,004 | 176 | 141 | 0 | 0 | 0 | 317 |
| Uadt_01 | 213 | 492,514 | 120 | 93 | 0 | 0 | 0 | 213 |

### 7.2 Incomplete Or Missing Studies

| Study | Status | Samples | Variants |
| --- | --- | ---: | ---: |
| Clrt_01 | bed-missing, bim-missing, fam-missing, summary-missing | - | - |
