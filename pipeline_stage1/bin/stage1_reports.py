#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
STAGE_NAME = "stage1"
TASK_SEQUENCE = ["summary"]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def count_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def write_tsv(rows: list[dict[str, object]], path: Path, columns: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")


def stage1_studies(analysis_root: Path, studies_arg: str) -> list[str]:
    if studies_arg != "all":
        return [study.strip() for study in studies_arg.split(",") if study.strip()]
    studies: list[str] = []
    for stage_dir in sorted(analysis_root.glob("*/stage1")):
        study = stage_dir.parent.name
        if (stage_dir / f"{study}.fam").exists() or (stage_dir / f"{study}.bim").exists():
            studies.append(study)
    return studies


def fallback_png(path: Path) -> None:
    data = (
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8"
        "/x8AAwMCAO+/p9sAAAAASUVORK5CYII="
    )
    path.write_bytes(base64.b64decode(data))


def fallback_pdf(path: Path, study: str) -> None:
    text = f"Stage 1 report counts for {study}"
    path.write_text(
        "%PDF-1.4\n"
        "1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n"
        "2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n"
        "3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 300 120] "
        "/Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >> endobj\n"
        f"4 0 obj << /Length {len(text) + 38} >> stream\n"
        f"BT /F1 12 Tf 20 70 Td ({text}) Tj ET\n"
        "endstream endobj\n"
        "5 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj\n"
        "xref\n0 6\n0000000000 65535 f \n"
        "trailer << /Size 6 /Root 1 0 R >>\nstartxref\n0\n%%EOF\n"
    )


def write_counts_figure(study: str, samples: int | None, variants: int | None, out_base: Path) -> None:
    ensure_dir(out_base.parent)
    png_path = out_base.with_suffix(".png")
    pdf_path = out_base.with_suffix(".pdf")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        labels = ["Samples", "Variants"]
        values = [samples or 0, variants or 0]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(labels, values, color=["#355c7d", "#6c5b7b"])
        ax.set_title(f"{study} stage 1 handoff")
        ax.set_ylabel("Count")
        for idx, value in enumerate(values):
            ax.text(idx, value, f"{value:,}", ha="center", va="bottom", fontsize=9)
        fig.tight_layout()
        fig.savefig(png_path)
        fig.savefig(pdf_path)
        plt.close(fig)
    except Exception:
        fallback_png(png_path)
        fallback_pdf(pdf_path, study)


def generate_study_report(analysis_root: Path, study: str) -> dict[str, object]:
    stage_dir = analysis_root / study / STAGE_NAME
    report_dir = stage_dir / "report"
    tables_dir = ensure_dir(report_dir / "tables")
    figures_dir = ensure_dir(report_dir / "figures")
    flags_dir = ensure_dir(report_dir / "flags")
    manifests_dir = ensure_dir(report_dir / "manifests")

    bed = stage_dir / f"{study}.bed"
    bim = stage_dir / f"{study}.bim"
    fam = stage_dir / f"{study}.fam"
    summary = stage_dir / "summary.txt"
    samples = count_lines(fam)
    variants = count_lines(bim)
    summary_lines = count_lines(summary)

    metrics_path = tables_dir / f"{study}_stage1_metrics.tsv"
    metrics = {
        "study": study,
        "samples": samples if samples is not None else "",
        "variants": variants if variants is not None else "",
        "has_bed": int(bed.exists()),
        "has_bim": int(bim.exists()),
        "has_fam": int(fam.exists()),
        "summary_lines": summary_lines if summary_lines is not None else "",
    }
    metric_columns = ["study", "samples", "variants", "has_bed", "has_bim", "has_fam", "summary_lines"]
    write_tsv([metrics], metrics_path, metric_columns)

    figure_base = figures_dir / f"{study}_stage1_counts"
    write_counts_figure(study, samples, variants, figure_base)

    flags = []
    for path in (bed, bim, fam):
        if not path.exists():
            flags.append({"study": study, "level": "ERROR", "message": f"Missing handoff file: {path.name}"})
    if summary_lines in (None, 0):
        flags.append({"study": study, "level": "WARN", "message": "Missing or empty stage summary."})
    write_tsv(flags, flags_dir / f"{study}_stage1_flags.tsv", ["study", "level", "message"])

    assets = [
        {"study": study, "asset_type": "table", "path": metrics_path},
        {"study": study, "asset_type": "figure", "path": figure_base.with_suffix(".png")},
        {"study": study, "asset_type": "figure", "path": figure_base.with_suffix(".pdf")},
    ]
    write_tsv(assets, manifests_dir / f"{study}_stage1_assets.tsv", ["study", "asset_type", "path"])
    write_tsv(assets, report_dir / "stage1_report_assets.tsv", ["study", "asset_type", "path"])
    return metrics


def generate_cohort_report(analysis_root: Path, rows: list[dict[str, object]]) -> None:
    report_dir = analysis_root / "cohort" / STAGE_NAME / "report"
    tables_dir = ensure_dir(report_dir / "tables")
    manifests_dir = ensure_dir(report_dir / "manifests")
    flags_dir = ensure_dir(report_dir / "flags")
    columns = ["study", "samples", "variants", "has_bed", "has_bim", "has_fam", "summary_lines"]
    summary_path = tables_dir / "stage1_summary.tsv"
    write_tsv(rows, summary_path, columns)
    assets = [{"study": "cohort", "asset_type": "table", "path": summary_path}]
    write_tsv(assets, manifests_dir / "stage1_assets.tsv", ["study", "asset_type", "path"])
    write_tsv([], flags_dir / "stage1_flags.tsv", ["study", "level", "message"])


def run_all_reports(analysis_root: Path, studies_arg: str, force: bool = False) -> None:
    del force
    rows = [generate_study_report(analysis_root, study) for study in stage1_studies(analysis_root, studies_arg)]
    generate_cohort_report(analysis_root, rows)


def run_task(
    task_id: str,
    analysis_root: Path,
    studies_arg: str,
    reference_dir: Path | None = None,
    force: bool = False,
) -> None:
    del task_id, reference_dir
    run_all_reports(analysis_root, studies_arg, force)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create read-only stage-1 report artifacts.")
    parser.add_argument("--analysis-root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--studies", default="all")
    parser.add_argument("--reference-dir", default=str(REPO_ROOT / "data" / "reference" / "1000G"))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main(task_id: str | None = None) -> None:
    args = parse_args()
    run_task(task_id or "summary", Path(args.analysis_root).resolve(), args.studies, Path(args.reference_dir), args.force)


if __name__ == "__main__":
    main()
