# pipeline_stage2

Stage 2 converts the stage-1 hg38 PLINK handoff into phased and imputed per-chromosome VCF outputs.

This stage consumes `analysis/<STUDY>/stage1/` outputs and writes the final imputed files to `analysis/<STUDY>/stage2/`. Temporary Nextflow work files are written under `pipeline_stage2/work/`.

## 1. Stage 2 Scope

Stage 2 starts from:

- `analysis/<STUDY>/stage1/<STUDY>.bed`
- `analysis/<STUDY>/stage1/<STUDY>.bim`
- `analysis/<STUDY>/stage1/<STUDY>.fam`

It also requires:

- the 1000 Genomes high-coverage GRCh38 reference panel under `data/reference/1000G/`
- the GRCh38 FASTA used for normalization
- the Eagle recombination map under `data/reference/eagle/genetic_map_hg38_withX.txt.gz`

Stage 2 ends with:

- `analysis/<STUDY>/stage2/<STUDY>_chr1_GxS.imputed.vcf.gz`
- `analysis/<STUDY>/stage2/<STUDY>_chr2_GxS.imputed.vcf.gz`
- `...`
- `analysis/<STUDY>/stage2/<STUDY>_chrX_GxS.imputed.vcf.gz` when chrX exists in stage 1

## 2. How To Run Stage 2

Full run:

```bash
sbatch src/005_stage2.sh
```

Single study:

```bash
sbatch src/005_stage2.sh --study Brea_01_Erneg
```


## 3. Ordered Execution Inside Stage 2

Stage 2 is defined in `pipeline_stage2/main.nf` and uses five active Nextflow modules:

- `reference_prep.nf`
- `target_prep.nf`
- `phasing.nf`
- `imputation.nf`
- `chrX_imputation.nf`

### 3.1 Step 1: Discover Stage-1 Inputs

The workflow scans:

- `analysis/*/stage1/*.{bed,bim,fam}`

It groups the three PLINK files into one logical stage-1 dataset per study and optionally filters that set using `--study`.

### 3.2 Step 2: Prepare The Reference Panel Per Chromosome

`PREP_REFERENCE` creates the imputation reference inputs chromosome by chromosome.

For autosomes, this module:

1. reads the 1000G phased reference VCF
2. filters to MAC `>= 10`
3. splits multiallelic sites
4. normalizes against the GRCh38 FASTA
5. writes a BCF plus index
6. compresses the BCF to Minimac4 `msav` format

For chrX, the pipeline uses a dedicated three-block strategy:

1. `PAR1`
2. `nonPAR`
3. `PAR2`

ChrX reference prep also keeps only SNPs and indels because the mixed SNV/INDEL/SV chrX panel caused `Minimac4 --compress-reference` failures in earlier testing.

### 3.3 Step 3: Export Stage-1 PLINK Files To Target VCF

`PREP_TARGET_VCF` converts each stage-1 PLINK dataset into one target VCF per chromosome.

For each chromosome `1..22` and `X`, the module:

1. checks whether that chromosome is actually present in the stage-1 `.bim`
2. skips chromosomes that are absent in stage 1
3. exports the chromosome with `plink --recode vcf bgz`
4. renames chromosomes to UCSC style (`chr1`, `chr2`, ..., `chrX`)
5. rewrites the VCF `ID` field as `CHROM:POS:REF:ALT`
6. splits multiallelic variants
7. writes and indexes the final target VCF

This step is where stage-1 study-specific choices matter most. For example, if stage 1 deliberately removed chrX, then stage 2 now skips chrX cleanly instead of treating that study as a task failure.

### 3.4 Step 4: Join Target VCFs To Their Matching Reference Chromosomes

The pipeline flattens the per-study target VCF outputs, extracts the chromosome from the filename, and joins each target chromosome to the reference object prepared in Step 2.

The joined channel is then split into:

- autosomes
- chrX

### 3.5 Step 5: Phase Autosomes With Eagle

`PHASE_AUTOSOMES` phases the autosomal target VCFs with Eagle.

For each study/chromosome pair, it receives:

- target VCF + index
- reference BCF + index
- the Eagle genetic map

The module then runs Eagle with:

- `--vcfRef`
- `--vcfTarget`
- `--allowRefAltSwap`
- `--geneticMapFile`

The phased output is:

- `<STUDY>_chr<CHR>_GxS.phased.vcf.gz`

### 3.6 Step 6: Impute Autosomes With Minimac4

`IMPUTE_AUTOSOMES` uses the phased autosomal VCF and the per-chromosome `msav` reference to run Minimac4.

For each autosome, the module:

1. loads the study’s phased target VCF
2. loads the corresponding chromosome `msav`
3. runs Minimac4
4. writes the imputed final VCF
5. indexes the final VCF with `bcftools index -t`

The published output is:

- `analysis/<STUDY>/stage2/<STUDY>_chr<CHR>_GxS.imputed.vcf.gz`

### 3.7 Step 7: Handle chrX Separately

`IMPUTE_CHRX` applies a dedicated chrX strategy.

The chrX pipeline:

1. splits the stage-2 target VCF into `PAR1`, `nonPAR`, and `PAR2`
2. intersects each block with the prepared reference positions
3. skips a block if there are no overlapping target/reference variants
4. phases each non-empty block with Eagle
5. imputes each non-empty block with Minimac4
6. concatenates the imputed chrX blocks that exist
7. resets final chrX IDs to `CHROM:POS:REF:ALT`

If a study has no chrX at all in stage 1, then chrX is absent from stage 2 by design and is not treated as an error.

