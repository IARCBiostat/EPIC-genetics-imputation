#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import html
import math
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_STAGES = ("stage1", "stage2", "stage3")
RENDER_STAGES = ("stage2", "stage3")   # stage1 omitted from master report sections
COPIED_STAGES = ("stage2", "stage3")

_R2_THRESHOLD = os.environ.get("STAGE3_MIN_R2", "0.3")
_MAF_THRESHOLD = os.environ.get("STAGE3_MAF", "0.01")
MAX_TABLE_ROWS = 20
FULL_ROW_COUNT_BYTE_LIMIT = 5_000_000

EXCLUDED_ANALYSIS_DIRS = {"cohort", "deprecated", "report", "reports"}

# Figures excluded: matched by endswith on lowercased filename (item 2, 3)
EXCLUDED_FIGURE_SUFFIXES = frozenset({'stage3_counts.png', 'stage1_counts.png'})

# Matches filenames like chr1.something.tsv or chrX.something.tsv
CHR_FILENAME_RE = re.compile(r'^(chr(?:\d+|X|Y))\.(.+)$', re.IGNORECASE)

# Per-chromosome files containing per-variant rows — too large to merge
LARGE_PER_CHROM_SUFFIXES = frozenset({'r2_maf.tsv', 'empirical_metrics.tsv'})

# Per-chromosome groups excluded from stage 2 (items 6, 15, 16 from prior session)
STAGE2_EXCLUDED_CHR_GROUP_SUFFIXES = frozenset({
    'imputation_metrics.tsv',
    'r2_maf.tsv',
    'empirical_metrics.tsv',
    'phase_quality.tsv',
    'phasing_metrics.tsv',   # item 6
})

# Per-chromosome groups excluded from stage 3 (item 10)
STAGE3_EXCLUDED_CHR_GROUP_SUFFIXES = frozenset({
    'variant_metrics.tsv',
})

# Tables that must show all rows without truncation (item 9)
FULL_DISPLAY_TABLES = frozenset({'stage3_variant_filters.tsv'})

# Human-readable stage names (items 5, 7)
STAGE_DISPLAY_NAMES = {
    "stage1": "Stage 1",
    "stage2": "Phasing And Imputation",
    "stage3": "QC And Finalisation",
}

DOSE0_NOTE = (
    "Empirical Dosage R2 and Dose0 are only populated for genotyped variants that "
    "underwent leave-one-out empirical validation; imputed variants show NA for these columns."
)

STAGE_EXPLANATIONS: dict[str, str] = {
    "stage2": (
        "<p>This stage takes the Stage&nbsp;1 quality-controlled and reference-harmonised "
        "genotype array data as input, phases the haplotypes into fully-resolved allele "
        "pairs, and then imputes ungenotyped variants using a population reference panel. "
        "All processing is performed independently per chromosome; results are subsequently "
        "merged and reported both per-chromosome and in aggregate. Detailed per-study reports "
        "are available in the stage-specific HTML files linked in the section below.</p>"

        "<p><strong>Input data.</strong> PLINK binary files (.bed/.bim/.fam) in genome build "
        "GRCh38 produced by Stage&nbsp;1, containing genotyped SNP array variants that have "
        "passed sample and variant quality control, strand alignment, and allele harmonisation "
        "to the GRCh38 reference sequence.</p>"

        "<p><strong>Phasing.</strong> Haplotype phasing is performed with "
        "<strong>EAGLE2</strong> using the GRCh38 sex-averaged genetic map "
        "(<code>genetic_map_hg38_withX.txt.gz</code>). Phasing conditions on the 1000 "
        "Genomes Project Phase&nbsp;3 reference haplotypes to maximise accuracy, "
        "particularly for low-frequency variants.</p>"

        "<p><strong>Imputation.</strong> Phased haplotypes are imputed with "
        "<strong>MINIMAC4</strong> against the <strong>1000 Genomes Project Phase&nbsp;3 "
        "reference panel</strong> (GRCh38 assembly; "
        "<code>GCA_000001405.15_GRCh38_no_alt_analysis_set</code>). Dosage genotypes and "
        "per-variant imputation quality scores are written to per-chromosome VCF files. The "
        "primary quality metric is <strong>R&sup2; (Rsq)</strong>, the estimated squared "
        "correlation between imputed and true genotype dosages. R&sup2; near 1.0 indicates "
        "high-confidence imputation; R&sup2; near 0 means dosage is driven almost entirely "
        "by allele frequency and carries little information about individual genotypes.</p>"

        "<p><strong>Empirical validation.</strong> For genotyped variants, a leave-one-out "
        "empirical validation is performed: each variant is masked and re-imputed from "
        "surrounding phased haplotypes, allowing the imputed dosage to be compared against "
        "the observed genotype. This yields an <em>empirical dosage R&sup2;</em> (distinct "
        "from the theoretical Rsq) and a <em>Dose0</em> metric (proportion of samples "
        "assigned a homozygous-reference dosage call). Both columns are NA for "
        "imputed-only variants.</p>"
    ),

    "stage3": (
        "<p>This stage applies post-imputation quality-control filters to the merged imputed "
        "dataset, annotates variants with dbSNP identifiers, conducts sample-level QC across "
        "four orthogonal criteria, and writes the final PLINK2 dataset ready for downstream "
        "analysis.</p>"

        "<p><strong>Variant filtering.</strong> Variants are retained only if they satisfy "
        "all of the following criteria:</p>"
        "<ul>"
        f"<li><strong>Imputation quality: R&sup2; &ge; {html.escape(_R2_THRESHOLD)}</strong> "
        f"&mdash; variants with lower imputation confidence are excluded regardless of allele "
        f"frequency (pipeline parameter <code>STAGE3_MIN_R2</code>).</li>"
        f"<li><strong>Minor allele frequency (MAF): &ge; {html.escape(_MAF_THRESHOLD)}</strong> "
        f"in the study population "
        f"(pipeline parameter <code>STAGE3_MAF</code>).</li>"
        "<li><strong>Hardy&ndash;Weinberg equilibrium (HWE):</strong> variants with a "
        "significant departure from HWE (exact test, applied in unrelated founders) are "
        "excluded where HWE filtering is enabled.</li>"
        "</ul>"
        "<p>After filtering, variant identifiers are annotated using a <strong>dbSNP</strong> "
        "VCF reference file: genotyped variants receive their rsID where available, imputed "
        "variants are matched by chromosomal position and alleles. Variants absent from dbSNP "
        "retain a <em>chr:pos:ref:alt</em> identifier.</p>"

        "<p><strong>Sex verification.</strong> Reported sex is compared against inferred sex "
        "derived from the X-chromosome inbreeding coefficient (F-statistic; PLINK2 "
        "<code>--check-sex</code>). Samples with a discordant reported and inferred sex are "
        "flagged and excluded from the final dataset.</p>"

        "<p><strong>Relatedness filtering.</strong> Pairwise kinship coefficients are "
        "estimated genome-wide using <strong>KING</strong> on a linkage-disequilibrium-pruned "
        "set of autosomal SNPs. For each pair of related individuals with kinship above the "
        "configured threshold (approximately second-degree relatives or closer), the sample "
        "with the lower genotyping call rate is removed. Removal counts are shown in the "
        "Stage&nbsp;3 Sample Filters table.</p>"

        "<p><strong>Heterozygosity filtering.</strong> Genome-wide heterozygosity is computed "
        "on a pruned set of common autosomal SNPs. Samples whose heterozygosity rate deviates "
        "by more than three standard deviations from the study mean are excluded as likely "
        "genotyping failures, sample contamination, or undetected duplicates.</p>"

        "<p><strong>Ancestry QC.</strong> Principal components analysis (PCA) is performed "
        "using PLINK2 on a pruned set of autosomal variants in a merged dataset of study "
        "samples and <strong>1000 Genomes Project Phase&nbsp;3 reference individuals</strong> "
        "spanning five super-populations (AFR, AMR, EAS, EUR, SAS). Samples whose PC "
        "coordinates identify them as outliers relative to the reference population clusters "
        "are excluded. The PCA plot is available in the Figures section below.</p>"

        "<p><strong>Output.</strong> The final dataset is written in PLINK2 binary format "
        "(.pgen/.pvar/.psam) to <code>{study}/stage3/final/</code>, containing only variants "
        "and samples that pass all QC filters. Note that the total unique samples removed is "
        "reported as the count of distinct individuals excluded across all four QC categories; "
        "this is typically less than the sum of per-category counts because some samples are "
        "independently flagged by more than one criterion.</p>"
    ),
}

TABLE_EXPLANATIONS = [
    ("variant_summary.tsv", "Variant Summary: number of overlapping genotyped variants used as the imputation basis, imputed target variants, and total variants by chromosome."),
    ("genotyped_maf_summary.tsv", "MAF Of Genotyped Variants: count of genotyped variants by MAF bin with empirical dosage R2 and Dose0 summaries where validation was run."),
    ("imputed_maf_summary.tsv", "MAF Of Imputed Variants: count of imputed variants by MAF bin with imputation R2 summaries."),
    ("empirical_validation_summary.tsv", "Empirical Validation Summary: leave-one-out empirical dosage R2, Dose0, and dosage-bias summaries for genotyped variants."),
    ("empirical_validation_by_maf.tsv", "Empirical Validation By MAF: empirical validation metrics aggregated into MAF bins."),
    ("r2_summary.tsv", (
        "Imputation R2 Summary: distribution of imputation R2 (Rsq, information score) across all "
        "imputed variants, summarised by chromosome and overall. R2 is the estimated squared "
        "correlation between imputed and true genotype dosages and is the primary quality measure "
        "for imputation: values near 1 indicate high-confidence imputation, while values below the "
        "configured threshold (default 0.3) are filtered in Stage 3. This table reports per-chromosome "
        "variant counts and R2 quantiles to characterise imputation quality across the genome."
    )),
    ("r2_by_maf_bins.tsv", "R2 By MAF Bins: imputation R2 and variant counts aggregated by MAF bin."),
    ("stage3_filter_steps.tsv", "Stage 3 Filter Steps: ordered post-imputation filtering components with input, output, and filtered counts."),
    ("stage3_variant_filters.tsv", "Stage 3 Variant Filters: per-chromosome R2/MAF and HWE filtering counts."),
    ("stage3_sample_filters.tsv", "Stage 3 Sample Filters: sex, relatedness, heterozygosity, ancestry, and unique sample-removal counts."),
    ("sample_review.tsv", "Stage 3 Sample Review: pre-final sample count and sample-removal category counts from the sample-review process."),
    ("variant_metrics.tsv", "Stage 3 Per-Chromosome Variant Metrics: post-imputation variant filtering and ID annotation metrics per chromosome."),
    ("imputation_metrics.tsv", "Chromosome Imputation Metrics: per-chromosome imputation output metrics staged during Stage 2."),
    ("r2_maf.tsv", "Chromosome R2 And MAF Metrics: per-variant R2 and MAF values. " + DOSE0_NOTE),
    ("empirical_metrics.tsv", "Chromosome Empirical Validation Metrics: per-variant leave-one-out validation and Dose0 metrics. " + DOSE0_NOTE),
    ("empirical_summary.tsv", "Chromosome Empirical Validation Summary: per-chromosome empirical validation summaries. " + DOSE0_NOTE),
    ("phasing_metrics.tsv", "Chromosome Phasing Metrics: per-chromosome phasing input/output counts and status metrics."),
    ("flags.tsv", "Flags: warnings or threshold flags raised by the reporting step."),
]

FIGURE_EXPLANATIONS = [
    ("af_concordance.png", "Allele-Frequency Concordance: reference allele frequency is plotted against study allele frequency to check concordance."),
    ("empirical_validation_by_maf.png", "Empirical Validation By MAF: empirical dosage R2 and Dose0 behaviour for genotyped variants across MAF bins."),
    ("r2_by_chromosome_violin.png", "Imputed SNP R2 Across Chromosomes: violin plot showing the distribution of imputation R2 by chromosome."),
    ("r2_by_maf.png", "Imputed SNP R2 By MAF: imputation R2 and variant counts across MAF bins for imputed variants."),
    ("ancestry_pca.png", "Stage 3 Ancestry Review: PCA scatter plot used to identify ancestry outliers."),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create per-study cross-stage master reports.")
    parser.add_argument("--analysis-root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--studies", default="all")
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def project_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return path.name


def resolve_study_id(analysis_root: Path, requested: str) -> str:
    exact = analysis_root / requested
    if exact.is_dir():
        return requested
    matches = [
        path.name
        for path in analysis_root.iterdir()
        if path.is_dir() and path.name.lower() == requested.lower()
    ]
    return matches[0] if len(matches) == 1 else requested


def study_ids(analysis_root: Path, studies_arg: str) -> list[str]:
    if studies_arg != "all":
        return sorted({resolve_study_id(analysis_root, study.strip()) for study in studies_arg.split(",") if study.strip()})

    studies: set[str] = set()
    for stage in REPORT_STAGES:
        for stage_dir in analysis_root.glob(f"*/{stage}"):
            if not stage_dir.is_dir():
                continue
            study = stage_dir.parent.name
            if study not in EXCLUDED_ANALYSIS_DIRS:
                studies.add(study)
    return sorted(studies, key=str.lower)


def report_dir(analysis_root: Path, study: str, stage: str) -> Path | None:
    path = analysis_root / study / stage / "report"
    return path if path.is_dir() else None


def clean_existing_copies(destination_dir: Path, study: str) -> None:
    if not destination_dir.exists():
        return
    prefix = f"{study}.".lower()
    for path in destination_dir.glob("*.html"):
        if path.name.lower().startswith(prefix):
            path.unlink()


def copy_stage_reports(analysis_root: Path, studies: list[str]) -> dict[tuple[str, str], list[Path]]:
    copied: dict[tuple[str, str], list[Path]] = {}
    report_root = ensure_dir(analysis_root / "report")
    for path in report_root.rglob(".DS_Store"):
        path.unlink(missing_ok=True)

    for stage in COPIED_STAGES:
        destination_dir = ensure_dir(report_root / stage)
        for study in studies:
            clean_existing_copies(destination_dir, study)
            source = report_dir(analysis_root, study, stage)
            if source is None:
                continue

            copied_paths: list[Path] = []
            html_reports = sorted(path for path in source.glob("*.html") if "master" not in path.name.lower())
            for html_report in html_reports:
                destination = destination_dir / f"{study}.{html_report.name}"
                destination.write_bytes(html_report.read_bytes())
                copied_paths.append(destination)
            copied[(study, stage)] = copied_paths

    return copied


def read_tsv(path: Path, limit: int | None = None) -> tuple[list[str], list[dict[str, str]], int | None, bool]:
    if not path.exists():
        return [], [], 0, False

    rows: list[dict[str, str]] = []
    total = 0
    truncated = False
    count_full_file = limit is None or path.stat().st_size <= FULL_ROW_COUNT_BYTE_LIMIT
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        for row in reader:
            total += 1
            if limit is None or len(rows) < limit:
                rows.append(row)
            elif limit is not None and not count_full_file:
                truncated = True
                break
            elif limit is not None:
                truncated = True
    return columns, rows, None if truncated and not count_full_file else total, truncated


def load_first_row(path: Path) -> dict[str, str]:
    _, rows, _, _ = read_tsv(path, limit=1)
    return rows[0] if rows else {}


def load_total_row(path: Path, label_column: str) -> dict[str, str]:
    _, rows, _, truncated = read_tsv(path, limit=None)
    if truncated:
        return {}
    for row in rows:
        if row.get(label_column, "").strip().lower() == "total":
            return row
    return rows[-1] if rows else {}


def count_non_header_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def fmt(value: object) -> str:
    if value is None or value == "":
        return "NA"
    text = str(value)
    try:
        number = float(text.replace(",", ""))
    except ValueError:
        return html.escape(text)
    if number.is_integer():
        return f"{int(number):,}"
    return f"{number:,.4f}"


def stage_display(stage: str) -> str:
    return STAGE_DISPLAY_NAMES.get(stage, stage.replace("stage", "Stage "))


def stripped_study_name(name: str, study: str) -> str:
    lowered = name.lower()
    study_lower = study.lower()
    for delimiter in ("_", "."):
        prefix = f"{study_lower}{delimiter}"
        if lowered.startswith(prefix):
            return name[len(study) + 1 :]
    return name


def artifact_key(path: Path, stage_report_dir: Path, study: str) -> str:
    relative = path.relative_to(stage_report_dir)
    name = stripped_study_name(relative.name, study)
    return str(relative.parent / name).lower()


def is_study_prefixed(path: Path, study: str) -> bool:
    name = path.name.lower()
    study_lower = study.lower()
    return name.startswith(f"{study_lower}_") or name.startswith(f"{study_lower}.")


def dedupe_artifacts(paths: list[Path], stage_report_dir: Path, study: str) -> list[Path]:
    grouped: dict[str, Path] = {}
    for path in paths:
        key = artifact_key(path, stage_report_dir, study)
        current = grouped.get(key)
        if current is None:
            grouped[key] = path
            continue
        current_is_prefixed = is_study_prefixed(current, study)
        new_is_prefixed = is_study_prefixed(path, study)
        if current_is_prefixed and not new_is_prefixed:
            grouped[key] = path
        elif current_is_prefixed == new_is_prefixed and len(str(path)) < len(str(current)):
            grouped[key] = path
    return sorted(grouped.values(), key=artifact_sort_key)


def artifact_sort_key(path: Path) -> tuple[int, str]:
    lower = path.as_posix().lower()
    priorities = [
        "stage3_filter_steps.tsv",
        "stage3_variant_filters.tsv",
        "stage3_sample_filters.tsv",
        "variant_summary.tsv",
        "genotyped_maf_summary.tsv",
        "imputed_maf_summary.tsv",
        "empirical_validation_summary.tsv",
        "ancestry_pca.png",
    ]
    for index, token in enumerate(priorities):
        if lower.endswith(token):
            return (index, lower)
    return (len(priorities), lower)


def human_title(path: Path, stage_report_dir: Path, study: str) -> str:
    relative = path.relative_to(stage_report_dir)
    name = stripped_study_name(relative.with_suffix("").as_posix(), study)
    return name.replace("/", " / ").replace("_", " ").replace(".", " ").title()


def explanation_for(path: Path, explanations: list[tuple[str, str]], default: str) -> str:
    lower = path.name.lower()
    relative_lower = path.as_posix().lower()
    for token, explanation in explanations:
        if lower.endswith(token) or token in relative_lower:
            return explanation
    return default


def encode_image(path: Path) -> str:
    if not path.exists():
        return ""
    return base64.b64encode(path.read_bytes()).decode("ascii")


def render_simple_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    header = "".join(f"<th>{html.escape(column.replace('_', ' ').title())}</th>" for column in columns)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{fmt(row.get(column, ''))}</td>" for column in columns) + "</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def render_tsv_table(path: Path, stage_report_dir: Path, study: str) -> str:
    # item 9: tables in FULL_DISPLAY_TABLES are shown without row limit
    limit = None if path.name.lower() in FULL_DISPLAY_TABLES else MAX_TABLE_ROWS
    columns, rows, total, truncated = read_tsv(path, limit=limit)
    if not columns:
        return ""

    header = "".join(f"<th>{html.escape(column.replace('_', ' ').title())}</th>" for column in columns)
    body_lines = []
    for row in rows:
        body_lines.append("<tr>" + "".join(f"<td>{fmt(row.get(column, ''))}</td>" for column in columns) + "</tr>")
    if truncated:
        body_lines.append(
            f"<tr><td colspan=\"{len(columns)}\">Showing first {len(rows):,} rows. Full file: {html.escape(project_relative(path))}</td></tr>"
        )

    title = human_title(path, stage_report_dir, study)
    explanation = explanation_for(
        path,
        TABLE_EXPLANATIONS,
        "Report table produced by the stage-specific reporting workflow.",
    )
    caption = (
        f"<p class=\"caption\">Table. {html.escape(explanation)} "
        f"Source: {html.escape(project_relative(path))}.</p>"
    )
    return f"<h4>{html.escape(title)}</h4>{caption}<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_lines)}</tbody></table>"


def render_figure(path: Path, stage_report_dir: Path, study: str) -> str:
    encoded = encode_image(path)
    if not encoded:
        return ""
    title = human_title(path, stage_report_dir, study)
    explanation = explanation_for(
        path,
        FIGURE_EXPLANATIONS,
        "PNG figure produced by the stage-specific reporting workflow.",
    )
    caption = f"Figure. {explanation} Source: {project_relative(path)}."
    return (
        "<figure>"
        f"<img src=\"data:image/png;base64,{encoded}\" alt=\"{html.escape(title)}\">"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


# ---------------------------------------------------------------------------
# Table filtering helpers
# ---------------------------------------------------------------------------

def is_asset_or_manifest_table(path: Path, stage_report_dir: Path) -> bool:
    relative = path.relative_to(stage_report_dir)
    if relative.parts and relative.parts[0].lower() == 'manifests':
        return True
    name_lower = path.name.lower()
    return (
        name_lower.endswith('_assets.tsv')
        or name_lower.endswith('_report_assets.tsv')
        or name_lower == 'assets.tsv'
    )


def is_flag_table(path: Path, stage_report_dir: Path) -> bool:
    relative = path.relative_to(stage_report_dir)
    if relative.parts and relative.parts[0].lower() == 'flags':
        return True
    name_lower = path.name.lower()
    return name_lower.endswith('_flags.tsv') or name_lower == 'flags.tsv'


def should_exclude_table(path: Path, stage_report_dir: Path, stage: str) -> bool:
    if is_asset_or_manifest_table(path, stage_report_dir):
        return True
    name_lower = path.name.lower()
    if name_lower.endswith('phasing_quality_summary.tsv'):
        return True
    # item 8: stage3_metrics rendered as separate data source, not as a table
    if stage == 'stage3' and name_lower.endswith('stage3_metrics.tsv'):
        return True
    if stage == 'stage2':
        for suffix in (
            'af_concordance.tsv',
            'af_concordance_summary.tsv',
            'empirical_validation_by_variant.tsv',
            'r2_by_variant.tsv',
        ):
            if name_lower.endswith(suffix):
                return True
    return False


def flag_source_label(path: Path, study: str) -> str:
    stem = path.stem.lower()
    study_lower = study.lower()
    for delim in ('_', '.'):
        if stem.startswith(f'{study_lower}{delim}'):
            stem = stem[len(study_lower) + 1:]
            break
    stem = re.sub(r'[_.]flags$', '', stem)
    if re.match(r'^stage\d+$', stem):
        stem = ''
    return stem.replace('_', ' ').replace('.', ' ').strip().title() or 'General'


def render_consolidated_flags(flag_paths: list[Path], stage_report_dir: Path, study: str) -> str:
    seen_cols: list[str] = []
    seen_col_set: set[str] = set()
    flagged_rows: list[dict[str, str]] = []

    for path in sorted(flag_paths):
        cols, rows, _, _ = read_tsv(path, limit=None)
        if not rows:
            continue
        label = flag_source_label(path, study)
        for col in cols:
            if col not in seen_col_set:
                seen_col_set.add(col)
                seen_cols.append(col)
        for row in rows:
            flagged_rows.append({'_check': label, **row})

    if not flagged_rows:
        return ''

    display_cols = ['_check'] + seen_cols
    header = ''.join(
        f'<th>{html.escape("Check" if c == "_check" else c.replace("_", " ").title())}</th>'
        for c in display_cols
    )
    body = ''.join(
        '<tr>' + ''.join(f'<td>{fmt(row.get(c, ""))}</td>' for c in display_cols) + '</tr>'
        for row in flagged_rows
    )
    return (
        '<h4>Flags</h4>'
        '<p class="caption">Table. Warnings and threshold flags raised during reporting. '
        'Only checks with active flags are shown.</p>'
        f'<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>'
    )


# ---------------------------------------------------------------------------
# Per-chromosome table consolidation
# ---------------------------------------------------------------------------

def chr_sort_key(chrom: str) -> tuple:
    c = chrom.lower()
    if c.startswith('chr'):
        c = c[3:]
    try:
        return (0, int(c), '')
    except ValueError:
        return (1, 0, c)


def canonical_chr_key(path: Path, study: str) -> tuple[str, str] | None:
    name = path.name
    name_lower = name.lower()
    study_lower = study.lower()
    for delim in ('_', '.'):
        prefix = f'{study_lower}{delim}'
        if name_lower.startswith(prefix):
            name = name[len(prefix):]
            break
    m = CHR_FILENAME_RE.match(name)
    if not m:
        return None
    chrom = m.group(1).lower()
    suffix = m.group(2)
    return suffix, chrom


def group_chr_tables(
    tables: list[Path],
    stage_report_dir: Path,
    study: str,
) -> tuple[list[tuple[str, list[Path], bool]], list[Path]]:
    buckets: dict[tuple[str, str], list[tuple[str, Path]]] = defaultdict(list)
    non_chr: list[Path] = []

    for path in tables:
        result = canonical_chr_key(path, study)
        if result:
            suffix, chrom = result
            parent_key = path.parent.relative_to(stage_report_dir).as_posix()
            buckets[(parent_key, suffix)].append((chrom, path))
        else:
            non_chr.append(path)

    groups: list[tuple[str, list[Path], bool]] = []
    for (parent_key, suffix), items in buckets.items():
        if len(items) < 2:
            non_chr.extend(p for _, p in items)
            continue
        sorted_paths = [p for _, p in sorted(items, key=lambda x: chr_sort_key(x[0]))]
        is_large = suffix in LARGE_PER_CHROM_SUFFIXES
        group_key = f'{parent_key}/{suffix}'
        groups.append((group_key, sorted_paths, is_large))

    return groups, non_chr


def chr_group_title(group_key: str) -> str:
    parts = group_key.split('/')
    if parts and parts[0].lower() in ('tables', '.'):
        parts = parts[1:]
    if parts:
        parts[-1] = parts[-1].replace('.tsv', '')
    title = ' / '.join(p.replace('_', ' ').title() for p in parts if p)
    return f'{title} (All Chromosomes)'


def render_chr_group(
    group_key: str,
    paths: list[Path],
    stage_report_dir: Path,
    study: str,
    is_large: bool,
) -> str:
    suffix = group_key.rsplit('/', 1)[-1]
    explanation = explanation_for(
        paths[0],
        TABLE_EXPLANATIONS,
        "Report table produced by the stage-specific reporting workflow.",
    )
    n_chrom = len(paths)
    title = chr_group_title(group_key)
    source_dir = html.escape(project_relative(paths[0].parent))

    if is_large:
        caption = (
            f"<p class=\"caption\">Table. {html.escape(explanation)} "
            f"Per-variant data across {n_chrom} chromosomes. "
            f"Full files available in the stage-specific report at {source_dir}.</p>"
        )
        return f"<h4>{html.escape(title)}</h4>{caption}"

    all_columns: list[str] = []
    all_rows: list[dict[str, str]] = []
    for path in paths:
        cols, rows, _, _ = read_tsv(path, limit=None)
        if cols and not all_columns:
            all_columns = cols
        all_rows.extend(rows)

    if not all_columns:
        return ""

    header = "".join(
        f"<th>{html.escape(c.replace('_', ' ').title())}</th>" for c in all_columns
    )
    body_lines = [
        "<tr>" + "".join(f"<td>{fmt(row.get(c, ''))}</td>" for c in all_columns) + "</tr>"
        for row in all_rows
    ]
    caption = (
        f"<p class=\"caption\">Table. {html.escape(explanation)} "
        f"Source directory: {source_dir}.</p>"
    )
    return (
        f"<h4>{html.escape(title)}</h4>{caption}"
        f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_lines)}</tbody></table>"
    )


# ---------------------------------------------------------------------------
# Sample flow SVG chart (item 4)
# ---------------------------------------------------------------------------

def _nice_ticks(lo: float, hi: float, n: int = 5) -> list[int]:
    """Compute n nicely-rounded tick values covering [lo, hi]."""
    span = hi - lo
    if span <= 0:
        return [int(lo)]
    raw = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    nice = raw / mag
    if nice <= 1.5:
        step = mag
    elif nice <= 3:
        step = 2 * mag
    elif nice <= 7:
        step = 5 * mag
    else:
        step = 10 * mag
    lo_t = math.floor(lo / step) * step
    ticks = []
    v = lo_t
    while v <= hi + step * 0.01:
        if v >= lo - step * 0.01:
            ticks.append(int(round(v)))
        v += step
    return ticks


def _svg_line_chart(points: list[tuple[str, int]], y_label: str = '', title: str = '') -> str:
    """Render an SVG line chart of (label, count) points."""
    if len(points) < 2:
        return ''

    W, H = 720, 300
    ML, MR, MT, MB = 80, 24, 44, 88
    cw = W - ML - MR
    ch = H - MT - MB

    counts = [p[1] for p in points]
    n = len(points)
    y_max = max(counts)
    y_min = min(counts)
    y_range = y_max - y_min or 1

    y_hi = y_max + y_range * 0.22   # headroom for value labels
    y_lo = max(0, y_min - y_range * 0.05)
    y_span = y_hi - y_lo

    def sx(i: int) -> float:
        return ML + i * cw / (n - 1)

    def sy(v: float) -> float:
        return MT + ch * (1.0 - (v - y_lo) / y_span)

    yticks = [t for t in _nice_ticks(y_lo, y_hi) if y_lo - 0.5 <= t <= y_hi + 0.5]

    lines: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'style="width:100%;max-width:{W}px;display:block;margin:0 auto;'
        f'font-family:Arial,sans-serif;">',
    ]

    if title:
        lines.append(
            f'<text x="{W / 2:.0f}" y="20" text-anchor="middle" '
            f'font-size="13" font-weight="bold" fill="#111827">{html.escape(title)}</text>'
        )

    # Grid and Y ticks
    for v in yticks:
        y = sy(v)
        lines.append(
            f'<line x1="{ML}" y1="{y:.1f}" x2="{ML + cw}" y2="{y:.1f}" '
            f'stroke="#e5e7eb" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{ML - 8}" y="{y:.1f}" text-anchor="end" '
            f'dominant-baseline="middle" font-size="11" fill="#6b7280">{v:,}</text>'
        )

    # Axes
    lines.append(
        f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT + ch}" stroke="#9ca3af" stroke-width="1.5"/>'
    )
    lines.append(
        f'<line x1="{ML}" y1="{MT + ch}" x2="{ML + cw}" y2="{MT + ch}" stroke="#9ca3af" stroke-width="1.5"/>'
    )

    # Y axis label
    if y_label:
        cx = 14
        cy = MT + ch // 2
        lines.append(
            f'<text transform="rotate(-90,{cx},{cy})" x="{cx}" y="{cy}" '
            f'text-anchor="middle" font-size="11" fill="#6b7280">{html.escape(y_label)}</text>'
        )

    # Line connecting points
    poly = ' '.join(f'{sx(i):.1f},{sy(c):.1f}' for i, c in enumerate(counts))
    lines.append(
        f'<polyline points="{poly}" fill="none" stroke="#3b82f6" '
        f'stroke-width="2.5" stroke-linejoin="round"/>'
    )

    # Dots, value labels, and X axis labels
    for i, (label, count) in enumerate(points):
        x = sx(i)
        y = sy(count)
        lines.append(
            f'<text x="{x:.1f}" y="{y - 10:.1f}" text-anchor="middle" '
            f'font-size="10" font-weight="bold" fill="#1d4ed8">{count:,}</text>'
        )
        lines.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="#3b82f6" stroke="#fff" stroke-width="2"/>'
        )
        lines.append(
            f'<text x="{x:.1f}" y="{MT + ch + 10}" text-anchor="end" '
            f'font-size="11" fill="#374151" '
            f'transform="rotate(-40,{x:.1f},{MT + ch + 10})">{html.escape(label)}</text>'
        )

    lines.append('</svg>')
    return '\n'.join(lines)


def generate_sample_flow_svg(analysis_root: Path, study: str, stage1_samples: int | None) -> str:
    """
    Build and return an SVG line chart of sample count through QC pipeline steps.
    Uses stage3_metrics.tsv for per-filter removal counts.
    Returns empty string if required data is not available.
    """
    if not stage1_samples:
        return ''

    stage3_rep = report_dir(analysis_root, study, 'stage3')
    if not stage3_rep:
        return ''

    metrics_path = stage3_rep / 'tables' / 'stage3_metrics.tsv'
    if not metrics_path.exists():
        metrics_path = stage3_rep / 'tables' / f'{study}_stage3_metrics.tsv'
    m = load_first_row(metrics_path)
    if not m:
        return ''

    def gi(col: str) -> int:
        try:
            return int(str(m.get(col, 0)).replace(',', ''))
        except (ValueError, TypeError):
            return 0

    n = stage1_samples
    final = gi('samples')
    sex_rm = gi('sex_mismatch')
    rel_rm = gi('related_removed')
    het_rm = gi('heterozygosity_outliers')

    # Build sequential points; the final count is the authoritative unique-deduplicated value
    points: list[tuple[str, int]] = [
        ('Stage 1', n),
        ('Post Sex Check', n - sex_rm),
        ('Post Relatedness', n - sex_rm - rel_rm),
        ('Post Het QC', n - sex_rm - rel_rm - het_rm),
        ('Final', final),
    ]

    return _svg_line_chart(
        points,
        y_label='Samples',
        title='Sample Count Through QC Pipeline',
    )


# ---------------------------------------------------------------------------
# Stage artifact discovery
# ---------------------------------------------------------------------------

def stage_artifacts(analysis_root: Path, study: str, stage: str) -> tuple[Path | None, list[Path], list[Path], list[Path]]:
    directory = report_dir(analysis_root, study, stage)
    if directory is None:
        return None, [], [], []

    html_reports = sorted(path for path in directory.glob("*.html") if "master" not in path.name.lower())
    figures = [path for path in directory.rglob("*.png") if path.name != ".DS_Store"]
    tables = [path for path in directory.rglob("*.tsv") if path.name != ".DS_Store"]
    figures = dedupe_artifacts(figures, directory, study)
    tables = dedupe_artifacts(tables, directory, study)
    return directory, html_reports, figures, tables


def overview_metrics(analysis_root: Path, study: str) -> dict[str, object]:
    # item 1: only final values are shown in the overview cards
    stage1_fam = analysis_root / study / "stage1" / f"{study}.fam"
    stage1_samples = count_non_header_lines(stage1_fam)

    stage3_rep = report_dir(analysis_root, study, "stage3")
    m: dict[str, str] = {}
    if stage3_rep:
        m = load_first_row(stage3_rep / "tables" / "stage3_metrics.tsv")
        if not m:
            m = load_first_row(stage3_rep / "tables" / f"{study}_stage3_metrics.tsv")

    final_samples_raw = m.get("samples", "")
    final_variants_raw = m.get("variants", "")

    try:
        final_s = int(str(final_samples_raw).replace(",", ""))
    except (ValueError, TypeError):
        final_s = None

    if stage1_samples and final_s is not None:
        total_removed: object = stage1_samples - final_s
    else:
        total_removed = m.get("total_removed", "")

    return {
        "stage1_samples": stage1_samples,   # kept for chart generation
        "final_samples": final_samples_raw,
        "final_variants": final_variants_raw,
        "total_samples_removed": total_removed,
    }


def render_stage_section(
    analysis_root: Path,
    copied: dict[tuple[str, str], list[Path]],
    study: str,
    stage: str,
) -> str:
    directory, html_reports, figures, tables = stage_artifacts(analysis_root, study, stage)
    if directory is None:
        return (
            f"<section><h2>{stage_display(stage)}</h2>"
            + STAGE_EXPLANATIONS.get(stage, "")
            + f"<p class=\"note\">No {stage_display(stage)} report directory was found for this study.</p>"
            "</section>"
        )

    # items 2, 3: exclude stage1_counts.png and stage3_counts.png
    filtered_figures = [
        p for p in figures
        if not any(p.name.lower().endswith(s) for s in EXCLUDED_FIGURE_SUFFIXES)
    ]
    figure_html = "\n".join(render_figure(path, directory, study) for path in filtered_figures)

    # Partition tables: flags consolidated separately; excluded tables dropped
    flag_tables: list[Path] = []
    normal_tables: list[Path] = []
    for path in tables:
        if is_flag_table(path, directory):
            flag_tables.append(path)
        elif not should_exclude_table(path, directory, stage):
            normal_tables.append(path)

    # Group per-chromosome tables, then apply stage-specific chr group exclusions
    chr_groups, non_chr_tables = group_chr_tables(normal_tables, directory, study)
    if stage == 'stage2':
        chr_groups = [
            (gk, paths, large) for gk, paths, large in chr_groups
            if gk.rsplit('/', 1)[-1] not in STAGE2_EXCLUDED_CHR_GROUP_SUFFIXES
        ]
    if stage == 'stage3':
        chr_groups = [
            (gk, paths, large) for gk, paths, large in chr_groups
            if gk.rsplit('/', 1)[-1] not in STAGE3_EXCLUDED_CHR_GROUP_SUFFIXES
        ]

    table_parts = [render_tsv_table(path, directory, study) for path in non_chr_tables]
    table_parts += [
        render_chr_group(gk, paths, directory, study, large)
        for gk, paths, large in chr_groups
    ]
    flags_html = render_consolidated_flags(flag_tables, directory, study)
    if flags_html:
        table_parts.append(flags_html)

    table_html = "\n".join(filter(None, table_parts))

    return f"""
<section>
<h2>{stage_display(stage)}</h2>
{STAGE_EXPLANATIONS.get(stage, '')}
<p class="note">Source report directory: <code>{html.escape(project_relative(directory))}</code>.</p>
<h3>Figures</h3>
{figure_html if figure_html else '<p class="note">No figures were found for this stage.</p>'}
<h3>Tables And Metrics</h3>
{table_html if table_html else '<p class="note">No tables were found for this stage.</p>'}
</section>
"""


def write_master_report(analysis_root: Path, copied: dict[tuple[str, str], list[Path]], study: str) -> Path:
    report_root = ensure_dir(analysis_root / "report")
    output_path = report_root / f"{study}.master-report.html"
    metrics = overview_metrics(analysis_root, study)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # item 4: sample flow SVG for the overview section
    sample_flow_svg = generate_sample_flow_svg(
        analysis_root, study, metrics["stage1_samples"]
    )

    stage_sections = "\n".join(
        render_stage_section(analysis_root, copied, study, stage) for stage in RENDER_STAGES
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Master Genetics Report: {html.escape(study)}</title>
<style>
body {{ margin: 0; background: #f3f4f6; color: #111827; font-family: Arial, sans-serif; line-height: 1.45; }}
.page {{ max-width: 1180px; margin: 0 auto; background: #fff; padding: 30px 36px 48px; }}
h1 {{ font-size: 30px; margin: 0 0 6px; }}
h2 {{ border-bottom: 2px solid #111827; font-size: 21px; margin-top: 32px; padding-bottom: 6px; }}
h3 {{ font-size: 17px; margin-top: 22px; }}
h4 {{ font-size: 14px; margin: 16px 0 4px; }}
.note, .caption, figcaption {{ color: #4b5563; font-size: 14px; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin: 18px 0 22px; }}
.metric {{ border: 1px solid #d1d5db; background: #fafafa; padding: 12px; }}
.metric span {{ display: block; color: #4b5563; font-size: 12px; text-transform: uppercase; }}
.metric strong {{ display: block; font-size: 22px; margin-top: 4px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 12px; margin: 8px 0 18px; }}
th, td {{ border: 1px solid #d1d5db; padding: 5px 7px; text-align: left; vertical-align: top; }}
th {{ background: #eef2f7; }}
figure {{ margin: 18px 0 24px; }}
img {{ display: block; width: min(100%, 900px); max-height: 620px; object-fit: contain; margin: 0 auto; border: 1px solid #d1d5db; }}
figcaption {{ max-width: 900px; margin: 7px auto 0; }}
ul {{ padding-left: 20px; }}
code {{ color: #374151; }}
footer {{ border-top: 1px solid #d1d5db; color: #6b7280; font-size: 12px; margin-top: 34px; padding-top: 12px; }}
@media (max-width: 760px) {{ .page {{ padding: 18px; }} img {{ width: 100%; max-height: none; }} }}
</style>
</head>
<body>
<main class="page">
<header>
<h1>Master Genetics Report: {html.escape(study)}</h1>
<p class="note">Generated: {generated}</p>
<p>This report collates the per-study pipeline outputs and QC metrics across the Phasing And Imputation (Stage&nbsp;2) and QC And Finalisation (Stage&nbsp;3) stages. Stage-specific HTML reports are copied into <code>analysis/report/stage2/</code> and <code>analysis/report/stage3/</code> and are cross-referenced below.</p>
</header>

<section>
<h2>Study Overview</h2>
<div class="metrics">
<div class="metric"><span>Final Samples</span><strong>{fmt(metrics["final_samples"])}</strong></div>
<div class="metric"><span>Final Variants</span><strong>{fmt(metrics["final_variants"])}</strong></div>
<div class="metric"><span>Total Samples Removed</span><strong>{fmt(metrics["total_samples_removed"])}</strong></div>
</div>
{sample_flow_svg}
<p class="caption" style="max-width:720px;margin:6px auto 0;">Figure. Sample count at each Stage 3 QC step. Points show the number of samples remaining after each filter. The final count reflects unique sample removal across all QC categories.</p>
</section>

{stage_sections}

<footer>Master report written to {html.escape(project_relative(output_path))}.</footer>
</main>
</body>
</html>
"""
    output_path.write_text(html_text)
    return output_path


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    if not analysis_root.is_dir():
        raise SystemExit(f"Analysis root does not exist: {analysis_root}")

    studies = study_ids(analysis_root, args.studies)
    if not studies:
        raise SystemExit("No studies were found for master report generation.")

    copied = copy_stage_reports(analysis_root, studies)
    for study in studies:
        output_path = write_master_report(analysis_root, copied, study)
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
