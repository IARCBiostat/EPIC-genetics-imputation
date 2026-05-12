#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import html
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_STAGES = ("stage1", "stage2", "stage3")
COPIED_STAGES = ("stage2", "stage3")
MAX_TABLE_ROWS = 20
FULL_ROW_COUNT_BYTE_LIMIT = 5_000_000

EXCLUDED_ANALYSIS_DIRS = {"cohort", "deprecated", "report", "reports"}

STAGE_EXPLANATIONS = {
    "stage1": (
        "Stage 1 standardises the raw genotype study, applies the study-specific "
        "pipeline filtering and harmonisation steps, and writes the authoritative "
        "Stage 1 PLINK handoff used by later stages."
    ),
    "stage2": (
        "Stage 2 phases the Stage 1 genotype data, imputes against the configured "
        "reference panel, and reports imputation performance. The detailed Stage 2 "
        "report describes phasing, imputation basis and target counts, R2 behaviour, "
        "allele-frequency concordance, and empirical validation where available."
    ),
    "stage3": (
        "Stage 3 performs post-imputation variant filtering and sample review. It "
        "applies configured R2, MAF, and HWE filters, reviews sex, relatedness, "
        "heterozygosity, and ancestry outliers, and writes the final PLINK2 handoff."
    ),
}

TABLE_EXPLANATIONS = [
    ("variant_summary.tsv", "Variant Summary: number of overlapping genotyped variants used as the imputation basis, imputed target variants, and total variants by chromosome."),
    ("genotyped_maf_summary.tsv", "MAF Of Genotyped Variants: count of genotyped variants by MAF bin with empirical dosage R2 and Dose0 summaries where validation was run."),
    ("imputed_maf_summary.tsv", "MAF Of Imputed Variants: count of imputed variants by MAF bin with imputation R2 summaries."),
    ("phasing_quality_summary.tsv", "Phasing Quality Summary: per-chromosome phasing metric availability and confidence summaries."),
    ("empirical_validation_summary.tsv", "Empirical Validation Summary: leave-one-out empirical dosage R2, Dose0, and dosage-bias summaries for genotyped variants."),
    ("empirical_validation_by_maf.tsv", "Empirical Validation By MAF: empirical validation metrics aggregated into MAF bins."),
    ("empirical_validation_by_variant.tsv", "Empirical Validation By Variant: per-variant leave-one-out empirical dosage R2, Dose0, dosage bias, and imputation R2 for genotyped validation variants."),
    ("af_concordance_summary.tsv", "Allele-Frequency Concordance Summary: study versus reference allele-frequency concordance metrics."),
    ("af_concordance.tsv", "Allele-Frequency Concordance: per-variant study and reference allele-frequency comparison values."),
    ("r2_summary.tsv", "R2 Summary: imputation R2 distribution metrics for imputed variants."),
    ("r2_by_maf_bins.tsv", "R2 By MAF Bins: imputation R2 and variant counts aggregated by MAF bin."),
    ("r2_by_variant.tsv", "R2 By Variant: per-variant imputation R2 values. Large files are previewed only in the master report."),
    ("stage3_metrics.tsv", "Stage 3 Metrics: final sample and variant counts plus aggregate Stage 3 filtering counts."),
    ("stage3_filter_steps.tsv", "Stage 3 Filter Steps: ordered post-imputation filtering components with input, output, and filtered counts."),
    ("stage3_variant_filters.tsv", "Stage 3 Variant Filters: per-chromosome R2/MAF and HWE filtering counts."),
    ("stage3_sample_filters.tsv", "Stage 3 Sample Filters: sex, relatedness, heterozygosity, ancestry, and unique sample-removal counts."),
    ("sample_review.tsv", "Stage 3 Sample Review: pre-final sample count and sample-removal category counts from the sample-review process."),
    ("variant_metrics.tsv", "Stage 3 Per-Chromosome Variant Metrics: post-imputation variant filtering and ID annotation metrics for one chromosome."),
    ("imputation_metrics.tsv", "Chromosome Imputation Metrics: per-chromosome imputation output metrics staged during Stage 2."),
    ("r2_maf.tsv", "Chromosome R2 And MAF Metrics: per-chromosome variant R2 and MAF values used by Stage 2 reports."),
    ("empirical_metrics.tsv", "Chromosome Empirical Validation Metrics: per-variant leave-one-out validation and Dose0 metrics."),
    ("empirical_summary.tsv", "Chromosome Empirical Validation Summary: per-chromosome empirical validation summaries."),
    ("phase_quality.tsv", "Chromosome Phase Quality: phase-quality values or availability flags for phased variants."),
    ("phasing_metrics.tsv", "Chromosome Phasing Metrics: per-chromosome phasing input/output counts and status metrics."),
    ("flags.tsv", "Flags: warnings or threshold flags raised by the reporting step."),
    ("assets.tsv", "Report Assets: manifest of report tables, figures, and HTML files produced for the stage."),
]

FIGURE_EXPLANATIONS = [
    ("af_concordance.png", "Allele-Frequency Concordance: reference allele frequency is plotted against study allele frequency to check concordance."),
    ("empirical_validation_by_maf.png", "Empirical Validation By MAF: empirical dosage R2 and Dose0 behaviour for genotyped variants across MAF bins."),
    ("r2_by_chromosome_violin.png", "Imputed SNP R2 Across Chromosomes: violin plot showing the distribution of imputation R2 by chromosome."),
    ("r2_by_maf.png", "Imputed SNP R2 By MAF: imputation R2 and variant counts across MAF bins for imputed variants."),
    ("stage3_counts.png", "Stage 3 Filtering Counts: variant and sample counts before and after Stage 3 filtering."),
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
        return str(path)


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
    return stage.replace("stage", "Stage ")


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
        "stage3_metrics.tsv",
        "stage3_variant_filters.tsv",
        "stage3_sample_filters.tsv",
        "variant_summary.tsv",
        "genotyped_maf_summary.tsv",
        "imputed_maf_summary.tsv",
        "phasing_quality_summary.tsv",
        "empirical_validation_summary.tsv",
        "af_concordance_summary.tsv",
        "stage3_counts.png",
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
    columns, rows, total, truncated = read_tsv(path, limit=MAX_TABLE_ROWS)
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

    total_text = f"{total:,}" if total is not None else f"at least {len(rows):,}"
    title = human_title(path, stage_report_dir, study)
    explanation = explanation_for(
        path,
        TABLE_EXPLANATIONS,
        "Report table produced by the stage-specific reporting workflow.",
    )
    caption = (
        f"<p class=\"caption\">Table. {html.escape(explanation)} "
        f"Rows: {total_text}. Source: {html.escape(project_relative(path))}.</p>"
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


def inventory_rows(analysis_root: Path, copied: dict[tuple[str, str], list[Path]], study: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for stage in REPORT_STAGES:
        directory, html_reports, figures, tables = stage_artifacts(analysis_root, study, stage)
        copied_reports = copied.get((study, stage), [])
        status = "present" if directory else "missing"
        if stage in COPIED_STAGES and directory and not html_reports:
            status = "missing stage HTML"
        rows.append(
            {
                "stage": stage_display(stage),
                "status": status,
                "source_report_dir": project_relative(directory) if directory else "",
                "copied_html_reports": len(copied_reports),
                "figures": len(figures),
                "tables": len(tables),
            }
        )
    return rows


def overview_metrics(analysis_root: Path, study: str) -> dict[str, object]:
    stage1_fam = analysis_root / study / "stage1" / f"{study}.fam"
    stage1_bim = analysis_root / study / "stage1" / f"{study}.bim"

    stage2_report = report_dir(analysis_root, study, "stage2")
    stage2_total = {}
    if stage2_report:
        stage2_total = load_total_row(stage2_report / "tables" / "variant_summary.tsv", "chromosome")

    stage3_report = report_dir(analysis_root, study, "stage3")
    stage3_metrics = {}
    if stage3_report:
        stage3_metrics = load_first_row(stage3_report / "tables" / "stage3_metrics.tsv")
        if not stage3_metrics:
            stage3_metrics = load_first_row(stage3_report / "tables" / f"{study}_stage3_metrics.tsv")

    return {
        "stage1_samples": count_non_header_lines(stage1_fam),
        "stage1_variants": count_non_header_lines(stage1_bim),
        "stage2_basis_variants": stage2_total.get("imputation_basis", ""),
        "stage2_imputed_variants": stage2_total.get("imputation_target", ""),
        "stage2_total_variants": stage2_total.get("total", ""),
        "stage3_samples": stage3_metrics.get("samples", ""),
        "stage3_variants": stage3_metrics.get("variants", ""),
        "stage3_samples_removed": stage3_metrics.get("total_removed", ""),
        "stage3_variants_removed": stage3_metrics.get("removed_r2_maf_variants", ""),
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
            f"<p>{html.escape(STAGE_EXPLANATIONS[stage])}</p>"
            f"<p class=\"note\">No {stage_display(stage)} report directory was found for this study.</p>"
            "</section>"
        )

    copied_reports = copied.get((study, stage), [])
    link_paths = copied_reports if copied_reports else html_reports
    links = ""
    if link_paths:
        link_items = "".join(
            f"<li><span>{html.escape(path.name)}</span><br><code>{html.escape(project_relative(path))}</code></li>"
            for path in link_paths
        )
        links = f"<h3>Study Report Files</h3><p class=\"caption\">Stage-specific HTML reports are copied centrally where applicable, while the source report directory remains the authority for tables and figures.</p><ul>{link_items}</ul>"

    figure_html = "\n".join(render_figure(path, directory, study) for path in figures)
    table_html = "\n".join(render_tsv_table(path, directory, study) for path in tables)

    return f"""
<section>
<h2>{stage_display(stage)}</h2>
<p>{html.escape(STAGE_EXPLANATIONS[stage])}</p>
<p class="note">Source report directory: <code>{html.escape(project_relative(directory))}</code>.</p>
{links}
<h3>Figures</h3>
{figure_html if figure_html else '<p class="note">No PNG figures were found for this stage.</p>'}
<h3>Tables, Metrics, Flags, And Manifests</h3>
{table_html if table_html else '<p class="note">No TSV tables were found for this stage.</p>'}
</section>
"""


def write_master_report(analysis_root: Path, copied: dict[tuple[str, str], list[Path]], study: str) -> Path:
    report_root = ensure_dir(analysis_root / "report")
    output_path = report_root / f"{study}.master-report.html"
    metrics = overview_metrics(analysis_root, study)
    inventory = inventory_rows(analysis_root, copied, study)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    stage_sections = "\n".join(render_stage_section(analysis_root, copied, study, stage) for stage in REPORT_STAGES)

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
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px; margin: 18px 0 22px; }}
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
<p>This report collates the per-study artifacts created across pipeline stages. Stage 2 and Stage 3 HTML reports are copied into <code>analysis/report/stage2/</code> and <code>analysis/report/stage3/</code> before the master report is written. Large TSV files are previewed, with the source path provided for the complete artifact.</p>
</header>

<section>
<h2>Study Overview</h2>
<div class="metrics">
<div class="metric"><span>Stage 1 Samples</span><strong>{fmt(metrics["stage1_samples"])}</strong></div>
<div class="metric"><span>Stage 1 Variants</span><strong>{fmt(metrics["stage1_variants"])}</strong></div>
<div class="metric"><span>Stage 2 Basis Variants</span><strong>{fmt(metrics["stage2_basis_variants"])}</strong></div>
<div class="metric"><span>Stage 2 Imputed Variants</span><strong>{fmt(metrics["stage2_imputed_variants"])}</strong></div>
<div class="metric"><span>Stage 2 Total Variants</span><strong>{fmt(metrics["stage2_total_variants"])}</strong></div>
<div class="metric"><span>Final Samples</span><strong>{fmt(metrics["stage3_samples"])}</strong></div>
<div class="metric"><span>Final Variants</span><strong>{fmt(metrics["stage3_variants"])}</strong></div>
<div class="metric"><span>Stage 3 Samples Removed</span><strong>{fmt(metrics["stage3_samples_removed"])}</strong></div>
<div class="metric"><span>Stage 3 R2/MAF Variant Removals</span><strong>{fmt(metrics["stage3_variants_removed"])}</strong></div>
</div>
</section>

<section>
<h2>Report Inventory</h2>
<p class="caption">This table checks whether each stage has a source report directory, whether copied HTML reports were created for Stage 2 and Stage 3, and how many figures and tables are included below.</p>
{render_simple_table(inventory, ["stage", "status", "source_report_dir", "copied_html_reports", "figures", "tables"])}
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
