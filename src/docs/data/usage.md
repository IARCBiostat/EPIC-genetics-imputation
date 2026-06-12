# Using Finalised Data

Each study's deliverable is written to `final/<STUDY>/` by stage 4.

## Directory Layout

```
final/<STUDY>/
├── <STUDY>.tar.gz        # per-chromosome PLINK2 files (pgen/pvar/psam)
├── report-master.html    # cross-stage master QC report
├── report-stage2.html    # stage 2 imputation report
├── report-stage3.html    # stage 3 sample and variant QC report
└── review/
    ├── related.exclude   # sample IDs flagged as related (KING ≥ 0.0884)
    ├── ancestry.exclude  # sample IDs flagged as ancestry outliers (PCA z-score ≥ 6.0)
    ├── hwe.exclude       # variant IDs failing HWE (p < 0.000005) on autosomes
    ├── het.het           # per-sample autosomal heterozygosity statistics
    ├── king.kin0         # pairwise kinship coefficients (KING format)
    ├── pca.eigenvec      # per-sample PCA scores (10 PCs)
    └── pca.eigenval      # PCA eigenvalues
```

## Extracting the Archive

```bash
cd final/<STUDY>
tar -xzf <STUDY>.tar.gz
```

The extracted directory contains per-chromosome PLINK2 files:

```
<STUDY>/
├── <STUDY>_chr1.pgen / .pvar
├── ...
├── <STUDY>_chr22.pgen / .pvar
├── <STUDY>_chrX.pgen / .pvar
└── <STUDY>.psam
```

Each chromosome is a self-contained PLINK2 dataset sharing the same `.psam` sample order.

## Loading in PLINK2

```bash
plink2 --pfile <STUDY>/<STUDY>_chr1 [...]
```

## Applying QC Exclusions

Heterozygosity outliers are removed from the PGEN files by default. Relatedness and ancestry exclusions are **off by default** and provided as `review/` files for analysts to apply as appropriate. HWE exclusions are also provided rather than hard-applied, as HWE filtering decisions are study- and analysis-specific.

### Sample exclusions

```bash
# Relatedness only
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --remove review/related.exclude \
  --make-pgen \
  --out <STUDY>_chr1_unrelated

# Relatedness + ancestry
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --remove <(cat review/related.exclude review/ancestry.exclude | sort -u) \
  --make-pgen \
  --out <STUDY>_chr1_filtered
```

### Variant exclusions

```bash
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --exclude review/hwe.exclude \
  --make-pgen \
  --out <STUDY>_chr1_hwe
```

### Loop across all chromosomes

```bash
STUDY="Glbd_01"
for chr in {1..22} X; do
  plink2 \
    --pfile ${STUDY}/${STUDY}_chr${chr} \
    --remove review/related.exclude \
    --exclude review/hwe.exclude \
    --make-pgen \
    --out ${STUDY}_qc/${STUDY}_chr${chr}
done
```

## Review Files Reference

| File | Contents | Use |
| --- | --- | --- |
| `related.exclude` | Sample IDs with KING kinship ≥ 0.0884 | `--remove` to drop related individuals |
| `ancestry.exclude` | Ancestry outliers on 10 PCs (z-score ≥ 6.0) | `--remove` to restrict to primary ancestry cluster |
| `hwe.exclude` | Variants failing HWE p < 0.000005 on autosomes | `--exclude` to remove HWE-failing variants |
| `het.het` | Per-sample autosomal heterozygosity (PLINK2 format) | Inspect outlier distribution; samples > 3 SD already removed from PGEN |
| `king.kin0` | Pairwise kinship coefficients (KING format) | Inspect relatedness graph; pairs > 0.0884 appear in `related.exclude` |
| `pca.eigenvec` | Per-sample scores on 10 PCs | Stratification correction in association analyses |
| `pca.eigenval` | Variance explained by each PC | Scree plots; deciding how many PCs to include as covariates |
