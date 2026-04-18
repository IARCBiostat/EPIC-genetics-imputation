# EPIC Genetics Pipeline

This repository contains the EPIC genetics imputation workflow. The workflow is composed of 3-stages:

1. `pipeline_stage1/`: study-specific raw [genotype](https://en.wikipedia.org/wiki/Genotype) preprocessing and liftover to hg38
2. `pipeline_stage2/`: [phasing](https://en.wikipedia.org/wiki/Haplotype_estimation) and [imputation](https://en.wikipedia.org/wiki/Imputation_(genetics)) against the [1000 Genomes Project](https://en.wikipedia.org/wiki/1000_Genomes_Project) high-coverage GRCh38 reference
3. `pipeline_stage3/`: post-imputation QC, [rsID](https://www.ncbi.nlm.nih.gov/snp/docs/RefSNP_about/) annotation, and conversion to PLINK2 format

The active Slurm submission wrappers are:

- `src/003_stage1.sh`
- `src/004_stage2.sh`
- `src/005_stage3.sh`

Stage-specific documentation lives in:

- [pipeline_stage1/README.md](pipeline_stage1/README.md)
- [pipeline_stage2/README.md](pipeline_stage2/README.md)
- [pipeline_stage3/README.md](pipeline_stage3/README.md)


## Stage 1: Study Harmonization and Genome-Build Standardization

- Directory: `pipeline_stage1/`
- Main purpose: Convert each raw study into a harmonized hg38 PLINK handoff
- Primary input:  Archive-style raw study PLINK data plus manifests and EPIC ID linkage files
- Primary output: `analysis/<STUDY>/stage1/<STUDY>.bed/.bim/.fam`
- Submission script: `src/003_stage1.sh` 

Stage 1 converts each raw study delivery into a common hg38 PLINK handoff. Each study is processed with its own bespoke script because the original deliveries differ in array platform, allele coding, chromosome naming, manifest structure, and identifier linkage rules. All downstream matching is coordinate-based and assumes one common genome build. All downstream tools expect each study to exist as one standardized PLINK handoff.

#### What this stage does:

- harmonizes sample identifiers to the EPIC identifier system
- applies study-specific preprocessing such as BIM sorting, AB-allele translation, sex-chromosome handling, and build-specific logic
- standardizes all studies onto GRCh38
- writes one consistent study-level PLINK dataset for the next stage

#### Data used:

- raw study genotype files from source
- study-specific chip manifests and metadata files
- study-specific EPIC linkage files
- the shared EPIC ID reference
- LiftOver resources used for build conversion to GRCh38

## Stage 2: Phasing and Imputation

- Directory: `pipeline_stage2/`
- Main purpose: Export stage-1 PLINK data to per-chromosome target VCFs, phase, and impute
- Primary input:  `analysis/<STUDY>/stage1/` plus 1000G GRCh38 reference and Eagle map
- Primary output: `analysis/<STUDY>/stage2/<STUDY>_chr*_GxS.imputed.vcf.gz`
- Submission script: `src/004_stage2.sh` 

#### What this stage does:
- performs phasing
- performs imputation

### Phasing

[Phasing](https://en.wikipedia.org/wiki/Haplotype_estimation) estimates which [alleles](https://en.wikipedia.org/wiki/Allele) sit together on the same physical [chromosome](https://en.wikipedia.org/wiki/Chromosome) copy. Standard [GWAS](https://en.wikipedia.org/wiki/Genome-wide_association_study) array data tell us which two alleles an individual has at a locus, but not which allele belongs to which parental [haplotype](https://en.wikipedia.org/wiki/Haplotype). Phasing resolves that ambiguity.

#### Why:

- imputation works best on phased haplotypes rather than unordered diploid genotypes
- phased target data align more naturally to phased reference haplotypes
- local haplotype structure carries more information than isolated variants

#### How:

- stage-1 PLINK files are exported to per-chromosome target VCFs
- Eagle phases each study chromosome against the matching public reference chromosome
- chrX is handled separately because `PAR1`, `nonPAR`, and `PAR2` require different treatment; see [here](https://github.com/statgen/Minimac4/issues/74)

#### Data used:

- stage-1 hg38 PLINK files, exported to chromosome-specific target VCFs
- the public 1000 Genomes high-coverage GRCh38 phased reference panel
	- the 1000 Genomes panel provides the external phased haplotypes that Eagle uses to infer study haplotypes
- the Eagle [genetic recombination map](https://en.wikipedia.org/wiki/Genetic_map)
	- the Eagle map supplies recombination information that helps decide where haplotypes are likely to switch

### Imputation

Imputation](https://en.wikipedia.org/wiki/Genotype_imputation) statistically infers [genotypes](https://en.wikipedia.org/wiki/Imputation_(genetics)) at variants that were not directly typed on the original array.

#### Why:

- different studies were genotyped on different arrays, so raw variant sets are not identical
- many informative variants are absent from one or more arrays
- imputation greatly increases genomic coverage and improves cross-study comparability

#### How:

- Minimac4 compares each phased study haplotype to the corresponding reference haplotypes
- when the study haplotypes closely match the reference haplotypes, untyped nearby alleles can be inferred
- the final output is an imputed per-chromosome VCF with dosage-style genotype information and quality metrics such as `R2`

#### Data used:

- phased study VCFs from the phasing step
- the same public 1000 Genomes GRCh38 reference panel, compressed to Minimac4 `msav` format
	- the 1000 Genomes panel is the public reference source that makes imputation possible
- the GRCh38 [reference genome](https://en.wikipedia.org/wiki/Reference_genome) FASTA used during normalization and reference preparation
	- the GRCh38 FASTA ensures that REF and ALT alleles are interpreted on the correct genome sequence

## Stage 3: Variant Annotation, QC, Conversion, and Sample QC

- Directory: `pipeline_stage3/`
- Main purpose: Annotate rsIDs, filter imputed variants, convert to PLINK2, and perform sample QC
- Primary input: `analysis/<STUDY>/stage2/` plus [dbSNP](https://en.wikipedia.org/wiki/DbSNP) GRCh38
- Primary output: `analysis/<STUDY>/stage3/<STUDY>.pgen/.pvar/.psam`
- Submission script: `src/005_stage3.sh`

#### What this stage does:

Stage 3 turns the chromosome-level imputed VCFs from stage 2 into one analysis-ready study-level PLINK2 dataset. It combines four linked operations into one workflow:

1. variant identifier standardization
2. variant-level QC
3. conversion to PLINK2 format
4. sample-level QC and final sample exclusion

#### Why:

- stage-2 outputs intentionally use coordinate-style IDs, but many downstream resources still expect rsIDs where available
- imputed data need an explicit post-imputation QC layer before they are treated as analysis-ready
- chromosome-level imputed VCFs are not the most practical final format for downstream study analyses
- sample-level QC is needed to identify technical problems, unexpected relatedness, sex mismatches, and ancestry outliers

#### How:

- each imputed chromosome is annotated against dbSNP
- if an exact GRCh38 match exists, the variant receives its rsID
- if no exact dbSNP match exists, the pipeline uses a fallback `CHR:POS:REF:ALT` identifier
- a mapping file is written so the identifier transition remains auditable
- the annotated variants are then filtered on imputation quality using `INFO/R2`
- the data are also filtered on [minor allele frequency](https://en.wikipedia.org/wiki/Minor_allele_frequency)
- [Hardy-Weinberg equilibrium](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle) filtering is applied by default on autosomes
- the filtered chromosome-level imputed VCFs are converted into PLINK2 `pgen/pvar/psam` files and merged into one study-level dataset
- the merged study-level dataset is then used for sample QC
- sample QC updates sex from the stage-1 `.fam`, estimates relatedness with KING, identifies [heterozygosity](https://en.wikipedia.org/wiki/Heterozygosity) outliers, performs sex checks, and runs [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)-based ancestry QC
- a final sample-removal list is built from the required QC exclusions and, optionally, ancestry outliers
- ancestry outliers are always identified
	- ancestry outliers are excluded only when the stage-3 exclusion flag is enabled

#### Data used:

- stage-2 imputed VCFs
- [dbSNP](https://en.wikipedia.org/wiki/DbSNP) GRCh38
	- dbSNP provides rsIDs where a known GRCh38 variant match exists
- the stage-1 `.fam` file as the authoritative recorded-sex source
	- the stage-1 `.fam` anchors sample sex to the original study handoff
- LD-pruned autosomal subsets derived from the stage-3 study data for relatedness and PCA-based QC
	- the LD-pruned autosomal subset provides a stable basis for relatedness, heterozygosity, and ancestry QC

## Output

Repository-level outputs are organized as:

- `analysis/<STUDY>/stage1/`
- `analysis/<STUDY>/stage2/`
- `analysis/<STUDY>/stage3/`
- `analysis/stage1-summary.md`
- `analysis/stage2-summary.md`
- `analysis/stage3-summary.md`

Temporary and workflow-specific files are kept inside the stage directories:

- `pipeline_stage1/work/`
- `pipeline_stage2/work/`
- `pipeline_stage3/work/`

## How To Run 

Prepare Tools And Reference Files

```bash
bash src/000_tools.sh
bash src/001_data-genetics.sh
bash src/002_data-reference.sh
```

Run Stages Sequentially

Full run:

```bash
sbatch src/003_stage1.sh
sbatch src/004_stage2.sh
sbatch src/005_stage3.sh
```

To exclude ancestry outliers during stage 3 finalization:

```bash
sbatch src/005_stage3.sh --exclude-ancestry-outliers
```


Test/Single study:

```bash
STAGE1_SCRIPTS=process_brea_01_erneg.py sbatch src/003_stage1.sh
sbatch src/004_stage2.sh --study Brea_01_Erneg
sbatch src/005_stage3.sh --study Brea_01_Erneg
```
To exclude ancestry outliers during stage 3 finalization:

```bash
sbatch src/005_stage3.sh --study Brea_01_Erneg --exclude-ancestry-outliers
```

## Methods And Thresholds Summary

### Stage 1

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| Study-specific preprocessing | Each study runs through one bespoke script to normalize file structure, SNP coding, and study-specific quirks before ID harmonization | No single global numeric threshold; logic is controlled by per-study settings for `BUILD`, `PRE`, `PART1`, `COMPLETION`, and `PART2_BUILD35_REL` | Methods include BIM sorting, AB-allele translation, name linkage, position reset, and optional early X/Y/XY/MT exclusion |
| Manifest comparison QC | `preProcessing.py` compares the raw study PLINK data to the array manifest and writes discrepancy files | No single global cut-off; exclusion sets are driven by manifest disagreement classes and archived study logic | Drives all downstream SNP exclusion and completion steps |
| Sample ID harmonization | Sample IDs are harmonized to EPIC IDs using the shared EPIC reference plus optional study-specific link files | `--update-sex`, `--remove`, `--update-ids`, `--make-pheno` with an empty phenotype file, `--indiv-sort n` | Stage 1 currently blanks phenotype intentionally after ID harmonization |
| SNP exclusion part 1 | Initial exclusion of manifest-disagreement SNPs | No single global numeric threshold; study-specific exclusion branches decide whether to remove negative-strand, chr/pos-mismatch, allele-mismatch, or unknown-strand SNPs | This is the first major study-specific branch point |
| Completion / repair step | Remaining variants have metadata repaired before the second exclusion pass | No single global numeric threshold; study-specific completion mode decides whether to update chromosome, position, strand, alleles, or rsID-based metadata | Converts the dataset into the best possible harmonized pre-final state |
| SNP exclusion part 2 | Final cleanup on the completed dataset | Always removes misplaced mitochondrial SNPs and duplicate classes 4, 2, 3, and 1; two studies also apply an extra build-35 exclusion list | Produces the final post-QC stage-1 prefix before liftover |
| Build harmonization | Final post-QC data are lifted to GRCh38 / hg38 | Build 36 -> `hg18ToHg38`; build 37 -> `hg19ToHg38`; final PLINK step uses `--split-x b38 no-fail` | Ensures all studies enter stage 2 on a common assembly |

### Stage 2

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| Target VCF export | Stage-1 PLINK data are exported chromosome by chromosome to target VCF | Chromosomes renamed to UCSC style; IDs rewritten to `CHROM:POS:REF:ALT`; multiallelics split | Chromosomes absent in stage 1 are skipped cleanly |
| Reference preparation | 1000 Genomes reference VCFs are normalized and converted to per-chromosome BCF / `msav` files | Reference filter includes `MAC >= 10` during prep; chrX keeps SNPs and indels only | chrX is handled as `PAR1`, `nonPAR`, and `PAR2` |
| Phasing | Eagle phases each study chromosome against the matching reference chromosome | `eagle_threads = 8`; Eagle uses the public recombination map | chrX is phased block by block rather than as one file |
| Imputation | Minimac4 imputes each phased chromosome against the matching `msav` reference panel | `min_r2 = 0.3`; `minimac_batch_size = 900`; `minimac_threads = 8` | Final stage-2 VCFs retain imputation INFO metrics including `R2` |
| chrX handling | chrX is processed separately from autosomes | Blocks: `PAR1`, `nonPAR`, `PAR2`; block is skipped if no overlap or too few target variants exist | Successful chrX blocks are concatenated; empty chrX output is written if none survive |
| Task runtime / retries | Internal stage-2 tasks are submitted to Slurm through Nextflow | Partition `low_p`; task walltime `10h`; retries only for exit codes `137`, `140`, `143`; `maxRetries = 2` | Intended to tolerate interruption-style failures without masking real data errors |

### Stage 3

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| rsID annotation | Each imputed chromosome is annotated against dbSNP GRCh38 | Exact dbSNP match -> keep rsID; otherwise fallback to `chr:pos:REF:ALT` | A per-chromosome mapping file is written for ID auditability |
| Variant QC: imputation quality | Variants are filtered on imputation quality | `min_r2 = 0.3` | Applied before PLINK2 conversion |
| Variant QC: allele frequency | Variants are filtered on minor allele frequency | `maf = 0.01` | Applied together with the `R2` filter |
| Variant QC: HWE | Hardy-Weinberg equilibrium filtering is applied by default on autosomes | `run_hwe = true`; `hwe_p = 0.000005`; `hwe_k = 0`; mode `midp keep-fewhet` | chrX is not HWE filtered |
| PLINK2 conversion | Filtered imputed VCFs are converted to chromosome-level `pgen/pvar/psam` and then merged | Uses dosage import from `HDS`; chrX uses `--split-par b38` and `--lax-chrx-import` | Produces one merged study-level PLINK2 dataset before sample QC |
| Sample QC: relatedness | Related or duplicate samples are identified with KING | `king_cutoff = 0.0884` | Related samples are always added to the removal list |
| Sample QC: heterozygosity | Heterozygosity outliers are identified from the LD-pruned autosomal dataset | `het_sd_threshold = 3.0` | Samples beyond the SD threshold are added to the removal list |
| Sample QC: ancestry identification | PCA-based ancestry outliers are identified for every study | `ancestry_pc_count = 10`; `ancestry_z_threshold = 6.0` | Identification is always performed |
| Sample QC: ancestry exclusion | Ancestry outliers may optionally be excluded from the final dataset | `exclude_ancestry_outliers = false` by default | Outliers are only removed when the exclusion flag is enabled |
| Sample QC: sex update and sex check | Recorded sex is updated from the stage-1 `.fam`, then genotype-derived sex is checked | No separate numeric threshold exposed in params | Sex mismatches are always added to the removal list |

## Study Overview Table

`XXX` indicates a value that will be populated from the stage-2 or stage-3 summaries after those summaries are regenerated. Raw-study counts are also currently left as `XXX` because there is not yet a consolidated raw-data summary file in the repository.

| Study | N raw samples | N raw variants | N stage 1 samples | N stage 1 variants | N stage 2 samples | N stage 2 variants | N stage 3 samples | N stage 3 variants |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Brea_01_Erneg` | XXX | XXX | 1,011 | 558,570 | XXX | XXX | XXX | XXX |
| `Brea_02` | XXX | XXX | 7,491 | 478,002 | XXX | XXX | XXX | XXX |
| `Clrt_01` | XXX | XXX | - | - | XXX | XXX | XXX | XXX |
| `Ecvd_01` | XXX | XXX | 9,426 | 531,696 | XXX | XXX | XXX | XXX |
| `Ecvd_02` | XXX | XXX | 8,920 | 208,221 | XXX | XXX | XXX | XXX |
| `Ecvd_03` | XXX | XXX | 8,587 | 395,482 | XXX | XXX | XXX | XXX |
| `Glbd_01` | XXX | XXX | 119 | 689,010 | XXX | XXX | XXX | XXX |
| `Inte_01` | XXX | XXX | 9,290 | 566,795 | XXX | XXX | XXX | XXX |
| `Inte_02` | XXX | XXX | 7,397 | 516,117 | XXX | XXX | XXX | XXX |
| `Inte_03` | XXX | XXX | 6,328 | 509,114 | XXX | XXX | XXX | XXX |
| `Kidn_01` | XXX | XXX | 356 | 580,897 | XXX | XXX | XXX | XXX |
| `Kidn_02` | XXX | XXX | 265 | 4,146,971 | XXX | XXX | XXX | XXX |
| `Lung_01` | XXX | XXX | 2,484 | 492,643 | XXX | XXX | XXX | XXX |
| `Lymp_01` | XXX | XXX | 480 | 732,277 | XXX | XXX | XXX | XXX |
| `Ovar_01` | XXX | XXX | 1,310 | 470,489 | XXX | XXX | XXX | XXX |
| `Panc_01` | XXX | XXX | 751 | 559,366 | XXX | XXX | XXX | XXX |
| `Panc_02` | XXX | XXX | 183 | 732,277 | XXX | XXX | XXX | XXX |
| `Pros_01` | XXX | XXX | 856 | 567,416 | XXX | XXX | XXX | XXX |
| `Pros_02` | XXX | XXX | 1,801 | 201,475 | XXX | XXX | XXX | XXX |
| `Pros_03` | XXX | XXX | 1,137 | 496,063 | XXX | XXX | XXX | XXX |
| `Pros_04` | XXX | XXX | 1,488 | 669,060 | XXX | XXX | XXX | XXX |
| `Stom_01` | XXX | XXX | 317 | 654,004 | XXX | XXX | XXX | XXX |
| `Uadt_01` | XXX | XXX | 213 | 492,514 | XXX | XXX | XXX | XXX |

## Metrics Table

This table mirrors the intended content of `analysis/stage3-summary.md`. For now it remains a scaffold and should be filled from the regenerated stage-3 summary once stage 3 has completed.

| Study | Variants In | Post R2/MAF | Post HWE | Final Variants | rsID | Fallback ID | Sex Mismatch | Related Removed | Het Outliers | Ancestry Identified | Ancestry Removed | Total Removed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Brea_01_Erneg` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Brea_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Clrt_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Ecvd_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Ecvd_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Ecvd_03` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Glbd_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Inte_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Inte_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Inte_03` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Kidn_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Kidn_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Lung_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Lymp_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Ovar_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Panc_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Panc_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Pros_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Pros_02` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Pros_03` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Pros_04` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Stom_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |
| `Uadt_01` | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX | XXX |

