# EPIC Genetics Pipeline

Multi-study GWAS imputation pipeline for the EPIC cancer cohort. The pipeline processes raw genotype array data through four stages — harmonisation to GRCh38, phasing and imputation against the 1000 Genomes reference, post-imputation QC and PLINK2 conversion, and final report generation — producing per-study analysis-ready datasets.

## Studies

| Study | N | Variants | Mean ER2 | Mean R2 | AF Pearson R |
| --- | ---: | ---: | ---: | ---: | ---: |
| Brea_01_Erneg | 987 | 11,593,992 | 0.9435 | 0.781 | 0.8729 |
| Brea_02 | 7,348 | 11,634,844 | 0.9328 | 0.685 | 0.8872 |
| Clrt_01 | 4,375 | 11,575,819 | 0.9513 | 0.761 | 0.8818 |
| Ecvd_01 | 9,238 | 11,617,557 | 0.9152 | 0.669 | 0.8873 |
| Ecvd_02 | 8,561 | 11,381,963 | 0.8923 | 0.496 | 0.9142 |
| Ecvd_03 | 8,479 | 9,644,865 | 0.9080 | 0.459 | 0.9149 |
| Glbd_01 | 114 | 11,617,976 | 0.9416 | 0.840 | 0.8518 |
| Inte_01 | 9,140 | 11,576,056 | 0.9436 | 0.733 | 0.8860 |
| Inte_02 | 7,244 | 11,648,231 | 0.9129 | 0.670 | 0.8885 |
| Inte_03 | 6,140 | 11,228,590 | 0.9197 | 0.685 | 0.8859 |
| Kidn_01 | 345 | 11,614,105 | 0.9420 | 0.814 | 0.8640 |
| Kidn_02 | 258 | 11,556,383 | 0.9778 | 0.909 | 0.8553 |
| Lung_01 | 2,414 | 11,655,446 | 0.9326 | 0.712 | 0.8846 |
| Lymp_01 | 457 | 11,659,973 | 0.9417 | 0.811 | 0.8711 |
| Neuro_01 | 4,830 | 11,852,050 | 0.9031 | 0.626 | 0.8941 |
| Ovar_01 | 1,282 | 11,688,077 | 0.9319 | 0.729 | 0.8826 |
| Panc_01 | 699 | 11,605,968 | 0.9427 | 0.794 | 0.8728 |
| Panc_02 | 174 | 11,585,552 | 0.9487 | 0.850 | 0.8549 |
| Pros_01 | 827 | 11,576,657 | 0.9461 | 0.791 | 0.8696 |
| Pros_03 | 1,115 | 11,658,676 | 0.9318 | 0.737 | 0.8802 |
| Pros_04 | 1,464 | 11,325,579 | 0.9339 | 0.756 | 0.8860 |
| Stom_01 | 308 | 11,632,113 | 0.9401 | 0.819 | 0.8646 |
| Uadt_01 | 206 | 11,652,795 | 0.9320 | 0.798 | 0.8630 |

*N and Variants reflect the Stage 3 final dataset after R²/MAF/HWE variant filtering and sample QC. Imputation metrics are from Stage 2.*

**Total unique participants across all 23 studies: 58,575** (76,005 total sample-study pairs; 11,416 participants appear in two or more studies).

## Sample Overlap

![Sample overlap across EPIC genetics studies](img/sample_overlap_upset.png)

