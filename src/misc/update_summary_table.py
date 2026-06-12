#!/usr/bin/env python3
"""
Update the README.md summary table and total-participants line from the latest pipeline outputs.

Reads (local):
  <summaries-dir>/stage3-summary.md  →  N (Stage3 Samples), Variants (Final Variants)
  <summaries-dir>/stage2-summary.md  →  Mean R2 (Mean Qual)
  <overlap-dir>/sample_overlap_summary.json  →  total unique, total pairs, n shared

Reads (scratch / analysis-root):
  <analysis-root>/<STUDY>/stage2/report/tables/empirical_validation_summary.tsv  →  Mean ER2
  <analysis-root>/<STUDY>/stage2/report/tables/af_concordance_summary.tsv        →  AF Pearson R

Writes:
  <readme>  updated in-place

Usage:
    python3 src/misc/update_summary_table.py \\
        --analysis-root /scratch/.../studies \\
        --summaries-dir analysis/final/summaries \\
        --overlap-dir report \\
        --readme README.md
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent.parent

TABLE_HEADER = "| Study | N | Variants | Mean ER2 | Mean R2 | AF Pearson R |"
DATA_PAGE_HEADER = "| Study | N | Variants | Mean ER2 | Report |"
TOTAL_PARTICIPANTS_RE = re.compile(
    r"\*\*Total unique participants across all \d+ studies: [\d,]+\*\*"
    r" \([\d,]+ total sample-study pairs; [\d,]+ participants appear in two or more studies\)\."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Update README.md summary table from latest pipeline outputs."
    )
    parser.add_argument(
        "--analysis-root",
        required=True,
        help="Path to the studies directory on scratch (e.g. ${SCRATCH_RUN}/studies)",
    )
    parser.add_argument(
        "--summaries-dir",
        default=str(REPO_ROOT / "analysis" / "final" / "summaries"),
        help="Path to local summaries directory containing stage2/3-summary.md",
    )
    parser.add_argument(
        "--overlap-dir",
        default=str(REPO_ROOT / "report"),
        help="Directory containing sample_overlap_summary.json (written by sample_overlap.py)",
    )
    parser.add_argument(
        "--readme",
        default=str(REPO_ROOT / "README.md"),
        help="Path to README.md to update in-place",
    )
    parser.add_argument(
        "--data-page",
        default=None,
        help=(
            "Optional path to the docs data page (src/docs/data/index.md). When "
            "given, its | Study | N | Variants | Mean ER2 | Report | table and "
            "total-participants line are regenerated from the same final-dataset source."
        ),
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_md_table(path: Path) -> list[dict[str, str]]:
    """Parse a markdown table into a list of row dicts keyed by header name."""
    if not path.exists():
        return []
    headers: list[str] = []
    rows: list[dict[str, str]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped.startswith("|"):
            continue
        cells = [c.strip() for c in stripped.strip("|").split("|")]
        if all(set(c) <= set("-: ") for c in cells):  # separator row
            continue
        if not headers:
            headers = cells
        else:
            rows.append(dict(zip(headers, cells)))
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fmt_int(value: str | int) -> str:
    try:
        return f"{int(str(value).replace(',', '')):,}"
    except (ValueError, TypeError):
        return str(value)


# ---------------------------------------------------------------------------
# Data readers
# ---------------------------------------------------------------------------

def read_stage3_stats(summaries_dir: Path) -> dict[str, dict[str, str]]:
    rows = parse_md_table(summaries_dir / "stage3-summary.md")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        study = row.get("Study", "").strip()
        if not study:
            continue
        result[study] = {
            "n": row.get("Stage3 Samples", "").replace(",", ""),
            "variants": row.get("Final Variants", "").replace(",", ""),
        }
    return result


def read_stage2_r2(summaries_dir: Path) -> dict[str, str]:
    rows = parse_md_table(summaries_dir / "stage2-summary.md")
    result: dict[str, str] = {}
    for row in rows:
        study = row.get("Study", "").strip()
        qual = row.get("Mean Qual", "").strip()
        if study and qual:
            result[study] = qual
    return result


def read_empirical_er2(study_dir: Path) -> str:
    """Compute weighted mean ER2 from per-chromosome empirical_validation_summary.tsv."""
    path = study_dir / "stage2" / "report" / "tables" / "empirical_validation_summary.tsv"
    rows = read_tsv(path)
    if not rows:
        return ""
    total_n = 0
    total_weighted = 0.0
    for row in rows:
        try:
            n = int(str(row.get("empirical_r2_variants", "0")).replace(",", ""))
            er2 = float(row.get("mean_empirical_dosage_r2", ""))
            if n > 0:
                total_n += n
                total_weighted += er2 * n
        except (ValueError, TypeError):
            continue
    if total_n == 0:
        return ""
    return f"{total_weighted / total_n:.4f}"


def read_pearson_r(study_dir: Path) -> str:
    path = study_dir / "stage2" / "report" / "tables" / "af_concordance_summary.tsv"
    rows = read_tsv(path)
    if not rows:
        return ""
    r = rows[0].get("pearson_r", "").strip()
    try:
        return f"{float(r):.4f}"
    except (ValueError, TypeError):
        return r


def read_overlap_summary(overlap_dir: Path) -> dict | None:
    path = overlap_dir / "sample_overlap_summary.json"
    if not path.exists():
        print(f"WARNING: sample_overlap_summary.json not found at {path}", file=sys.stderr)
        return None
    with path.open() as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# README update
# ---------------------------------------------------------------------------

def build_table_rows(
    studies: list[str],
    stage3: dict[str, dict[str, str]],
    stage2_r2: dict[str, str],
    analysis_root: Path,
) -> list[str]:
    """README / home-page table: | Study | N | Variants | Mean ER2 | Mean R2 | AF Pearson R |"""
    rows = []
    missing_er2: list[str] = []
    missing_pearson: list[str] = []
    for study in studies:
        s3 = stage3.get(study, {})
        n = fmt_int(s3.get("n", ""))
        variants = fmt_int(s3.get("variants", ""))
        r2 = stage2_r2.get(study, "")
        study_dir = analysis_root / study
        er2 = read_empirical_er2(study_dir)
        pearson = read_pearson_r(study_dir)
        if not er2:
            missing_er2.append(study)
        if not pearson:
            missing_pearson.append(study)
        rows.append(f"| {study} | {n} | {variants} | {er2} | {r2} | {pearson} |")
    if missing_er2:
        print(f"WARNING: Mean ER2 missing for: {', '.join(missing_er2)}", file=sys.stderr)
    if missing_pearson:
        print(f"WARNING: AF Pearson R missing for: {', '.join(missing_pearson)}", file=sys.stderr)
    return rows


def build_data_page_rows(
    studies: list[str],
    stage3: dict[str, dict[str, str]],
    analysis_root: Path,
) -> list[str]:
    """Data-page table: | Study | N | Variants | Mean ER2 | Report |

    N and Variants come from the same Stage 3 final-dataset source as the README,
    so every N reported on the site is the final `.psam` sample count.
    """
    rows = []
    for study in studies:
        s3 = stage3.get(study, {})
        n = fmt_int(s3.get("n", ""))
        variants = fmt_int(s3.get("variants", ""))
        er2 = read_empirical_er2(analysis_root / study)
        report = f"[report](reports/{study}.master-report.html)"
        rows.append(f"| {study} | {n} | {variants} | {er2} | {report} |")
    return rows


def _build_total_line(overlap: dict) -> str:
    return (
        f"**Total unique participants across all {overlap.get('n_studies', 0)} studies: "
        f"{overlap.get('total_unique', 0):,}**"
        f" ({overlap.get('total_pairs', 0):,} total sample-study pairs;"
        f" {overlap.get('n_shared', 0):,} participants appear in two or more studies)."
    )


def _replace_table_and_total(
    path: Path,
    header: str,
    table_rows: list[str],
    overlap: dict | None,
) -> None:
    """Replace the markdown table under `header` and the total-participants line in `path`."""
    lines = path.read_text().splitlines(keepends=True)
    new_lines: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if header in line:
            new_lines.append(line)
            i += 1
            # Preserve separator row
            if i < len(lines) and lines[i].strip().startswith("|") and set(lines[i].replace("|", "").replace(" ", "").replace(":", "").replace("\n", "")) <= set("-"):
                new_lines.append(lines[i])
                i += 1
            # Skip old data rows
            while i < len(lines) and lines[i].strip().startswith("|"):
                i += 1
            # Insert new rows
            for row in table_rows:
                new_lines.append(row + "\n")
        else:
            new_lines.append(line)
            i += 1

    text = "".join(new_lines)
    if overlap:
        text = TOTAL_PARTICIPANTS_RE.sub(_build_total_line(overlap), text)
    path.write_text(text)


def update_readme(
    readme_path: Path,
    table_rows: list[str],
    overlap: dict | None,
) -> None:
    _replace_table_and_total(readme_path, TABLE_HEADER, table_rows, overlap)


def update_data_page(
    data_page_path: Path,
    table_rows: list[str],
    overlap: dict | None,
) -> None:
    _replace_table_and_total(data_page_path, DATA_PAGE_HEADER, table_rows, overlap)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root)
    summaries_dir = Path(args.summaries_dir)
    overlap_dir = Path(args.overlap_dir)
    readme_path = Path(args.readme)

    for label, path in [
        ("analysis root", analysis_root),
        ("summaries dir", summaries_dir),
        ("README", readme_path),
    ]:
        if not path.exists():
            raise SystemExit(f"ERROR: {label} not found: {path}")

    stage3 = read_stage3_stats(summaries_dir)
    if not stage3:
        raise SystemExit("ERROR: no studies found in stage3-summary.md")

    stage2_r2 = read_stage2_r2(summaries_dir)
    overlap = read_overlap_summary(overlap_dir)
    studies = sorted(stage3.keys())

    table_rows = build_table_rows(studies, stage3, stage2_r2, analysis_root)
    update_readme(readme_path, table_rows, overlap)
    print(f"Updated {readme_path} with {len(studies)} studies.")

    if args.data_page:
        data_page_path = Path(args.data_page)
        if not data_page_path.exists():
            raise SystemExit(f"ERROR: data page not found: {data_page_path}")
        data_page_rows = build_data_page_rows(studies, stage3, analysis_root)
        update_data_page(data_page_path, data_page_rows, overlap)
        print(f"Updated {data_page_path} with {len(studies)} studies.")

    if overlap:
        print(
            f"Total unique participants: {overlap['total_unique']:,}"
            f" across {overlap['n_studies']} studies"
            f" ({overlap['total_pairs']:,} sample-study pairs)."
        )


if __name__ == "__main__":
    main()
