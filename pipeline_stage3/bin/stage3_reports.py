#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import html
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
STAGE_NAME = "stage3"
TASK_SEQUENCE = ["summary"]
EXPECTED_CHROMS = [str(chrom) for chrom in range(1, 23)] + ["X"]
CHROM_RE = re.compile(r"(?:^|_)chr([0-9]+|X)(?:[._]|$)", re.IGNORECASE)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Stage 3 post-imputation report artifacts.")
    parser.add_argument("--analysis-root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--studies", default="all")
    parser.add_argument("--reference-dir", default=str(REPO_ROOT / "data" / "reference" / "1000G"))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_params() -> dict[str, str]:
    params_path = REPO_ROOT / "pipeline_stage3" / "params.yaml"
    params: dict[str, str] = {}
    if not params_path.exists():
        return params
    with params_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or ":" not in stripped:
                continue
            key, value = stripped.split(":", 1)
            params[key.strip()] = value.strip().strip('"').strip("'")
    return params


def project_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path)


def count_non_header_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_first_tsv_row(path: Path) -> dict[str, str]:
    rows = read_tsv_rows(path)
    return rows[0] if rows else {}


def int_value(row: dict[str, object], key: str, default: int = 0) -> int:
    value = row.get(key, default)
    if value is None or value == "":
        return default
    try:
        return int(float(str(value).replace(",", "")))
    except ValueError:
        return default


def fmt_int(value: int | None) -> str:
    if value is None:
        return ""
    return f"{value:,}"


def fmt_cell(value: object) -> str:
    if value is None or value == "":
        return ""
    text = str(value)
    try:
        number = int(float(text.replace(",", "")))
    except ValueError:
        return html.escape(text)
    return f"{number:,}"


def write_tsv(rows: list[dict[str, object]], path: Path, columns: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def chrom_from_path(path: Path) -> str:
    match = CHROM_RE.search(path.name)
    return match.group(1).upper() if match else ""


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    chrom = chrom.replace("chr", "").replace("CHR", "")
    if chrom.isdigit():
        return (int(chrom), "")
    if chrom.upper() == "X":
        return (23, "")
    if chrom.upper() == "Y":
        return (24, "")
    return (99, chrom)


def study_ids(analysis_root: Path, studies_arg: str) -> list[str]:
    if studies_arg != "all":
        return sorted({study.strip() for study in studies_arg.split(",") if study.strip()})

    excluded = {"cohort", "deprecated", "report", "reports"}
    studies: set[str] = set()
    for stage_dir in analysis_root.glob("*/stage3"):
        if not stage_dir.is_dir():
            continue
        study = stage_dir.parent.name
        if study not in excluded:
            studies.add(study)
    return sorted(studies)


def cleanup_stale_report_outputs(report_dir: Path) -> None:
    if not report_dir.exists():
        return
    (report_dir / ".DS_Store").unlink(missing_ok=True)
    for path in report_dir.rglob(".DS_Store"):
        path.unlink(missing_ok=True)
    for path in report_dir.rglob("*.pdf"):
        path.unlink(missing_ok=True)
    for pattern in ("*master_report*.html", "*master-report*.html"):
        for path in report_dir.glob(pattern):
            path.unlink(missing_ok=True)


def fallback_png(path: Path) -> None:
    data = (
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8"
        "/x8AAwMCAO+/p9sAAAAASUVORK5CYII="
    )
    ensure_dir(path.parent)
    path.write_bytes(base64.b64decode(data))


def archive_path(stage_dir: Path, study: str) -> Path:
    return stage_dir / f"{study}.stage3.tar.gz"


def summarize_variant_metrics(tables_dir: Path) -> tuple[list[dict[str, object]], dict[str, int]]:
    rows: list[dict[str, object]] = []
    totals = {
        "input_variants": 0,
        "rsid_variants": 0,
        "fallback_variants": 0,
        "duplicate_rsid_fallbacks": 0,
        "post_r2_maf_variants": 0,
        "removed_r2_maf_variants": 0,
        "hwe_exclude_count": 0,
        "hwe_chromosomes": 0,
        "hwe_no_controls_chromosomes": 0,
    }

    for path in sorted(tables_dir.glob("*.variant_metrics.tsv"), key=lambda item: chrom_sort_key(chrom_from_path(item))):
        source = read_first_tsv_row(path)
        if not source:
            continue
        chrom = (source.get("chr") or chrom_from_path(path)).replace("chr", "")
        input_variants = int_value(source, "input_variants")
        rsid_variants = int_value(source, "rsid_variants")
        fallback_variants = int_value(source, "fallback_variants")
        duplicate_rsid_fallbacks = int_value(source, "duplicate_rsid_fallbacks")
        post_r2_maf = int_value(source, "post_r2_maf_variants")
        hwe_applied = int_value(source, "hwe_applied")
        hwe_no_controls = int_value(source, "hwe_no_controls")
        hwe_exclude = int_value(source, "hwe_exclude_count")
        removed_r2_maf = max(input_variants - post_r2_maf, 0)

        row = {
            "chromosome": chrom,
            "input_variants": input_variants,
            "rsid_variants": rsid_variants,
            "fallback_variants": fallback_variants,
            "duplicate_rsid_fallbacks": duplicate_rsid_fallbacks,
            "post_r2_maf_variants": post_r2_maf,
            "removed_r2_maf_variants": removed_r2_maf,
            "hwe_exclude_count": hwe_exclude,
            "hwe_applied": hwe_applied,
            "hwe_no_controls": hwe_no_controls,
        }
        rows.append(row)

        for key in totals:
            if key in ("hwe_chromosomes", "hwe_no_controls_chromosomes"):
                continue
            totals[key] += int_value(row, key)
        totals["hwe_chromosomes"] += 1 if hwe_applied else 0
        totals["hwe_no_controls_chromosomes"] += 1 if hwe_no_controls else 0

    return rows, totals


def build_filter_steps(
    variant_totals: dict[str, int],
    sample_row: dict[str, str],
    final_samples: int | None,
    params: dict[str, str],
) -> list[dict[str, object]]:
    pre_final_samples = int_value(sample_row, "pre_final_samples")
    total_removed = int_value(sample_row, "total_removed")
    min_r2 = params.get("min_r2", "0.3")
    maf = params.get("maf", "0.01")
    hwe_p = params.get("hwe_p", "0.000005")

    return [
        {
            "component": "Input handoff",
            "measurement": "variants",
            "input_count": variant_totals["input_variants"],
            "output_count": variant_totals["input_variants"],
            "filtered_count": 0,
            "details": "Imputed variants entering Stage 3 post-imputation filtering.",
        },
        {
            "component": "Variant ID annotation",
            "measurement": "variants",
            "input_count": variant_totals["input_variants"],
            "output_count": variant_totals["input_variants"],
            "filtered_count": 0,
            "details": "dbSNP rsIDs are retained where available; fallback chr:pos:ref:alt IDs are used otherwise.",
        },
        {
            "component": "Imputation quality and MAF filter",
            "measurement": "variants",
            "input_count": variant_totals["input_variants"],
            "output_count": variant_totals["post_r2_maf_variants"],
            "filtered_count": variant_totals["removed_r2_maf_variants"],
            "details": f"Configured filter: INFO/R2 >= {min_r2} and INFO/MAF >= {maf}.",
        },
        {
            "component": "Hardy-Weinberg exclusion list",
            "measurement": "variants",
            "input_count": variant_totals["post_r2_maf_variants"],
            "output_count": variant_totals["post_r2_maf_variants"],
            "filtered_count": variant_totals["hwe_exclude_count"],
            "details": f"HWE assessed on controls only (autosomes, p <= {hwe_p}; chromosome X excluded). Failing variants are written to hwe.exclude for downstream use — no hard filter applied.",
        },
        {
            "component": "Sample review",
            "measurement": "samples",
            "input_count": pre_final_samples,
            "output_count": final_samples if final_samples is not None else "",
            "filtered_count": total_removed,
            "details": "Union of relatedness, heterozygosity, and configured ancestry exclusions.",
        },
    ]


def build_sample_filter_rows(sample_row: dict[str, str]) -> list[dict[str, object]]:
    rows = [
        ("Relatedness identified", "related_identified", "Samples identified by KING relatedness cutoff (always computed)."),
        ("Relatedness removed", "related_removed", "Samples removed due to relatedness (only when --related true is set)."),
        ("Heterozygosity", "heterozygosity_outliers", "Samples outside the configured heterozygosity threshold."),
        ("Ancestry identified", "ancestry_outliers_identified", "Samples outside the configured PCA ancestry threshold."),
        ("Ancestry removed", "ancestry_outliers_removed", "Ancestry outliers added to the removal list when configured."),
        ("Total unique removed", "total_removed", "Unique sample IDs removed before the final Stage 3 handoff."),
    ]
    return [
        {"filter": label, "sample_count": int_value(sample_row, key), "details": details}
        for label, key, details in rows
    ]


def write_counts_figure(
    study: str,
    variant_totals: dict[str, int],
    pre_final_samples: int,
    final_samples: int | None,
    out_path: Path,
) -> None:
    ensure_dir(out_path.parent)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
        variant_labels = ["Input", "Post R2/MAF"]
        variant_values = [
            variant_totals["input_variants"],
            variant_totals["post_r2_maf_variants"],
        ]
        sample_labels = ["Pre-review", "Final"]
        sample_values = [pre_final_samples, final_samples or 0]

        axes[0].bar(variant_labels, variant_values, color=["#4b5563", "#059669"])
        axes[0].set_title("Variant filtering")
        axes[0].set_ylabel("Variants")
        axes[0].tick_params(axis="x", labelrotation=25)

        axes[1].bar(sample_labels, sample_values, color=["#4b5563", "#059669"])
        axes[1].set_title("Sample filtering")
        axes[1].set_ylabel("Samples")

        for ax in axes:
            for idx, value in enumerate([patch.get_height() for patch in ax.patches]):
                ax.text(idx, value, f"{int(value):,}", ha="center", va="bottom", fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.suptitle(f"{study} Stage 3 filtering counts", fontsize=12)
        fig.tight_layout()
        fig.savefig(out_path, dpi=160)
        plt.close(fig)
    except Exception:
        fallback_png(out_path)


def write_ancestry_figure(stage_dir: Path, report_dir: Path, study: str, out_path: Path) -> None:
    ensure_dir(out_path.parent)
    eigenvec_path = stage_dir / "sample_review" / f"{study}_pca.eigenvec"
    outliers_path = report_dir / "flags" / f"{study}.ancestry_outliers.id"
    legacy_png = report_dir / "figures" / f"{study}_ancestry_pca.png"

    if not eigenvec_path.exists() and legacy_png.exists() and legacy_png != out_path:
        shutil.copyfile(legacy_png, out_path)
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        if not eigenvec_path.exists():
            fallback_png(out_path)
            return

        df = pd.read_csv(eigenvec_path, sep=r"\s+")
        outliers = set()
        if outliers_path.exists():
            with outliers_path.open() as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) >= 2:
                        outliers.add(fields[1])

        df["stage3_status"] = df["IID"].map(lambda value: "Filtered" if value in outliers else "Retained")
        colors = {"Retained": "#2563eb", "Filtered": "#dc2626"}

        fig, ax = plt.subplots(figsize=(6.5, 5.2))
        for label, group in df.groupby("stage3_status"):
            ax.scatter(group["PC1"], group["PC2"], s=20, alpha=0.7, label=label, color=colors[label])
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("Stage 3 ancestry review")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(out_path, dpi=160)
        plt.close(fig)
    except Exception:
        fallback_png(out_path)


def render_table(rows: list[dict[str, object]], columns: list[str], max_rows: int | None = None) -> str:
    visible_rows = rows[:max_rows] if max_rows else rows
    header = "".join(f"<th>{html.escape(column.replace('_', ' ').title())}</th>" for column in columns)
    body_lines = []
    for row in visible_rows:
        body_lines.append("<tr>" + "".join(f"<td>{fmt_cell(row.get(column, ''))}</td>" for column in columns) + "</tr>")
    if max_rows and len(rows) > max_rows:
        body_lines.append(
            f"<tr><td colspan=\"{len(columns)}\">Showing {max_rows:,} of {len(rows):,} rows.</td></tr>"
        )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_lines)}</tbody></table>"


def encode_image(path: Path) -> str:
    if not path.exists():
        return ""
    return base64.b64encode(path.read_bytes()).decode("ascii")


def figure_html(path: Path, alt: str, caption: str) -> str:
    encoded = encode_image(path)
    if not encoded:
        return ""
    return (
        "<figure>"
        f"<img src=\"data:image/png;base64,{encoded}\" alt=\"{html.escape(alt)}\">"
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


def write_html_report(
    study: str,
    report_dir: Path,
    metrics_row: dict[str, object],
    filter_steps: list[dict[str, object]],
    variant_rows: list[dict[str, object]],
    sample_filter_rows: list[dict[str, object]],
    figure_paths: list[tuple[Path, str, str]],
) -> Path:
    output_path = report_dir / "report-stage3.html"
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    figures = "\n".join(figure_html(path, alt, caption) for path, alt, caption in figure_paths)

    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Stage 3 Report: {html.escape(study)}</title>
<style>
body {{ margin: 0; background: #f4f5f7; color: #17202a; font-family: Arial, sans-serif; line-height: 1.45; }}
.page {{ max-width: 1120px; margin: 0 auto; background: #fff; padding: 28px 34px 44px; }}
h1 {{ font-size: 28px; margin: 0 0 6px; }}
h2 {{ border-bottom: 2px solid #17202a; font-size: 20px; margin-top: 30px; padding-bottom: 6px; }}
h3 {{ font-size: 16px; margin: 18px 0 6px; }}
p {{ margin: 8px 0 12px; }}
.note, figcaption {{ color: #4b5563; font-size: 14px; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px; margin: 18px 0; }}
.metric {{ border: 1px solid #d8dee4; background: #fafafa; padding: 12px; }}
.metric span {{ display: block; color: #4b5563; font-size: 12px; text-transform: uppercase; }}
.metric strong {{ display: block; font-size: 22px; margin-top: 4px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin: 8px 0 16px; }}
th, td {{ border: 1px solid #d8dee4; padding: 6px 8px; text-align: left; vertical-align: top; }}
th {{ background: #eef2f7; }}
figure {{ margin: 18px 0 24px; }}
img {{ display: block; width: min(100%, 860px); max-height: 620px; object-fit: contain; margin: 0 auto; border: 1px solid #d8dee4; }}
figcaption {{ max-width: 860px; margin: 7px auto 0; }}
footer {{ border-top: 1px solid #d8dee4; color: #6b7280; font-size: 12px; margin-top: 34px; padding-top: 12px; }}
@media (max-width: 760px) {{ .page {{ padding: 18px; }} img {{ width: 100%; max-height: none; }} }}
</style>
</head>
<body>
<main class="page">
<header>
<h1>Stage 3 Report: {html.escape(study)}</h1>
<p class="note">Generated: {generated}</p>
<p>Stage 3 performs post-imputation filtering and sample review. It starts from the imputed VCFs, assigns final variant IDs, applies configured variant filters, reviews sample-level exclusions, and writes the final PLINK2 handoff files.</p>
</header>

<section>
<h2>Key Metrics</h2>
<div class="metrics">
<div class="metric"><span>Final Samples</span><strong>{fmt_cell(metrics_row.get("samples"))}</strong></div>
<div class="metric"><span>Final Variants</span><strong>{fmt_cell(metrics_row.get("variants"))}</strong></div>
<div class="metric"><span>R2/MAF Filtered Variants</span><strong>{fmt_cell(metrics_row.get("removed_r2_maf_variants"))}</strong></div>
<div class="metric"><span>HWE Exclusion List</span><strong>{fmt_cell(metrics_row.get("hwe_exclude_count"))}</strong></div>
<div class="metric"><span>Total Samples Removed</span><strong>{fmt_cell(metrics_row.get("total_removed"))}</strong></div>
</div>
</section>

<section>
<h2>Filtering Steps</h2>
<p class="note">Table 1. Stage 3 components in execution order. Counts show how many variants or samples entered each component, how many remained after it, and how many were filtered by that component.</p>
{render_table(filter_steps, ["component", "measurement", "input_count", "output_count", "filtered_count", "details"])}
</section>

<section>
<h2>Variant Filtering By Chromosome</h2>
<p class="note">Table 2. Per-chromosome post-imputation variant counts. The R2/MAF columns describe the imputation-quality and allele-frequency filter; hwe_exclude_count is the number of variants written to the per-study hwe.exclude file (no hard filter is applied).</p>
{render_table(variant_rows, ["chromosome", "input_variants", "rsid_variants", "fallback_variants", "duplicate_rsid_fallbacks", "post_r2_maf_variants", "removed_r2_maf_variants", "hwe_exclude_count", "hwe_applied", "hwe_no_controls"])}
</section>

<section>
<h2>Sample Review</h2>
<p class="note">Table 3. Sample-level Stage 3 review. Category counts can overlap; the total removed row is the unique set of samples removed before final handoff.</p>
{render_table(sample_filter_rows, ["filter", "sample_count", "details"])}
</section>

<section>
<h2>Figures</h2>
{figures}
</section>

<footer>Report artifacts are stored under {html.escape(project_relative(report_dir))}.</footer>
</main>
</body>
</html>
"""
    output_path.write_text(html_text)
    return output_path


def generate_study_report(analysis_root: Path, study: str) -> dict[str, object]:
    params = read_params()
    stage_dir = analysis_root / study / STAGE_NAME
    report_dir = ensure_dir(stage_dir / "report")
    tables_dir = ensure_dir(report_dir / "tables")
    figures_dir = ensure_dir(report_dir / "figures")
    flags_dir = ensure_dir(report_dir / "flags")
    manifests_dir = ensure_dir(report_dir / "manifests")
    cleanup_stale_report_outputs(report_dir)

    variant_rows, variant_totals = summarize_variant_metrics(tables_dir)
    sample_review_path = tables_dir / f"{study}.sample_review.tsv"
    sample_row = read_first_tsv_row(sample_review_path)

    pre_final_samples = int_value(sample_row, "pre_final_samples")
    total_removed = int_value(sample_row, "total_removed")
    samples = max(pre_final_samples - total_removed, 0) if pre_final_samples else None
    variants = variant_totals["post_r2_maf_variants"] or None

    metrics_row = {
        "study": study,
        "samples": samples if samples is not None else "",
        "variants": variants if variants is not None else "",
        "input_variants": variant_totals["input_variants"],
        "post_r2_maf_variants": variant_totals["post_r2_maf_variants"],
        "removed_r2_maf_variants": variant_totals["removed_r2_maf_variants"],
        "hwe_exclude_count": variant_totals["hwe_exclude_count"],
        "rsid_variants": variant_totals["rsid_variants"],
        "fallback_variants": variant_totals["fallback_variants"],
        "duplicate_rsid_fallbacks": variant_totals["duplicate_rsid_fallbacks"],
        "hwe_chromosomes": variant_totals["hwe_chromosomes"],
        "related_identified": sample_row.get("related_identified", ""),
        "related_removed": sample_row.get("related_removed", ""),
        "heterozygosity_outliers": sample_row.get("heterozygosity_outliers", ""),
        "ancestry_outliers_identified": sample_row.get("ancestry_outliers_identified", ""),
        "ancestry_outliers_removed": sample_row.get("ancestry_outliers_removed", ""),
        "total_removed": sample_row.get("total_removed", ""),
    }

    metrics_columns = list(metrics_row)
    write_tsv([metrics_row], tables_dir / "stage3_metrics.tsv", metrics_columns)
    write_tsv([metrics_row], tables_dir / f"{study}_stage3_metrics.tsv", metrics_columns)

    filter_steps = build_filter_steps(variant_totals, sample_row, samples, params)
    sample_filter_rows = build_sample_filter_rows(sample_row)
    variant_columns = [
        "chromosome",
        "input_variants",
        "rsid_variants",
        "fallback_variants",
        "duplicate_rsid_fallbacks",
        "post_r2_maf_variants",
        "removed_r2_maf_variants",
        "hwe_exclude_count",
        "hwe_applied",
        "hwe_no_controls",
    ]
    write_tsv(filter_steps, tables_dir / "stage3_filter_steps.tsv", ["component", "measurement", "input_count", "output_count", "filtered_count", "details"])
    write_tsv(variant_rows, tables_dir / "stage3_variant_filters.tsv", variant_columns)
    write_tsv(sample_filter_rows, tables_dir / "stage3_sample_filters.tsv", ["filter", "sample_count", "details"])

    counts_figure = figures_dir / "stage3_counts.png"
    ancestry_figure = figures_dir / "ancestry_pca.png"
    write_counts_figure(study, variant_totals, pre_final_samples, samples, counts_figure)
    write_ancestry_figure(stage_dir, report_dir, study, ancestry_figure)

    flags: list[dict[str, object]] = []
    if not variant_rows:
        flags.append({"study": study, "level": "WARN", "message": "No per-chromosome variant metrics were found."})
    if not sample_row:
        flags.append({"study": study, "level": "WARN", "message": "No sample review metrics were found."})
    study_archive = archive_path(stage_dir, study)
    if not study_archive.exists():
        flags.append({"study": study, "level": "WARN", "message": f"Stage 3 archive not found: {project_relative(study_archive)}"})
    write_tsv(flags, flags_dir / "stage3_flags.tsv", ["study", "level", "message"])
    write_tsv(flags, flags_dir / f"{study}_stage3_flags.tsv", ["study", "level", "message"])

    html_path = write_html_report(
        study,
        report_dir,
        metrics_row,
        filter_steps,
        variant_rows,
        sample_filter_rows,
        [
            (counts_figure, "Stage 3 filtering counts", "Figure 1. Stage 3 variant and sample counts before and after post-imputation filtering."),
            (ancestry_figure, "Stage 3 ancestry review", "Figure 2. Principal component scatter plot used for Stage 3 ancestry outlier review. Filtered ancestry outliers are highlighted when eigenvector data are available."),
        ],
    )

    assets = [
        {"study": study, "asset_type": "html", "path": project_relative(html_path)},
        {"study": study, "asset_type": "table", "path": project_relative(tables_dir / "stage3_metrics.tsv")},
        {"study": study, "asset_type": "table", "path": project_relative(tables_dir / "stage3_filter_steps.tsv")},
        {"study": study, "asset_type": "table", "path": project_relative(tables_dir / "stage3_variant_filters.tsv")},
        {"study": study, "asset_type": "table", "path": project_relative(tables_dir / "stage3_sample_filters.tsv")},
        {"study": study, "asset_type": "figure", "path": project_relative(counts_figure)},
        {"study": study, "asset_type": "figure", "path": project_relative(ancestry_figure)},
    ]
    write_tsv(assets, manifests_dir / "stage3_assets.tsv", ["study", "asset_type", "path"])
    write_tsv(assets, manifests_dir / f"{study}_stage3_assets.tsv", ["study", "asset_type", "path"])
    write_tsv(assets, report_dir / "stage3_report_assets.tsv", ["study", "asset_type", "path"])

    return metrics_row


def run_all_reports(analysis_root: Path, studies_arg: str, force: bool = False) -> None:
    del force
    for study in study_ids(analysis_root, studies_arg):
        generate_study_report(analysis_root, study)


def run_task(
    task_id: str,
    analysis_root: Path,
    studies_arg: str,
    reference_dir: Path | None = None,
    force: bool = False,
) -> None:
    del task_id, reference_dir
    run_all_reports(analysis_root, studies_arg, force)


def main(task_id: str | None = None) -> None:
    args = parse_args()
    run_task(task_id or "summary", Path(args.analysis_root).resolve(), args.studies, Path(args.reference_dir), args.force)


if __name__ == "__main__":
    main()
