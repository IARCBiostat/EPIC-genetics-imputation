#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


STUDY_ID_RE = re.compile(r"^STUDY_ID = '(.*)'$", re.MULTILINE)


@dataclass
class StudySummary:
    study_id: str
    bed_exists: bool
    bim_exists: bool
    fam_exists: bool
    summary_exists: bool
    sample_count: int | None
    sample_ids: frozenset[tuple[str, str]]
    variant_count: int | None
    sex_1: int
    sex_2: int
    sex_other: int
    pheno_1: int
    pheno_2: int
    pheno_other: int

    @property
    def complete(self) -> bool:
        return self.bed_exists and self.bim_exists and self.fam_exists and self.summary_exists

    @property
    def file_status(self) -> str:
        parts = []
        parts.append("bed" if self.bed_exists else "bed-missing")
        parts.append("bim" if self.bim_exists else "bim-missing")
        parts.append("fam" if self.fam_exists else "fam-missing")
        parts.append("summary" if self.summary_exists else "summary-missing")
        return ", ".join(parts)


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Write a markdown summary of stage-1 outputs.")
    parser.add_argument(
        "--analysis-root",
        default=str(repo_root / "analysis"),
        help="Root directory containing analysis/<STUDY>/stage1 outputs.",
    )
    parser.add_argument(
        "--scripts-dir",
        default=str(repo_root / "pipeline_stage1" / "scripts"),
        help="Directory containing process_*.py study scripts.",
    )
    parser.add_argument(
        "--output",
        default=str(repo_root / "analysis" / "stage1-summary.md"),
        help="Path to the markdown summary to write.",
    )
    return parser.parse_args()


def read_study_ids(scripts_dir: Path) -> list[str]:
    study_ids: list[str] = []
    for script_path in sorted(scripts_dir.glob("process_*.py")):
        text = script_path.read_text()
        match = STUDY_ID_RE.search(text)
        if match:
            study_ids.append(match.group(1))
    return study_ids


def count_lines(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for _ in handle)


def summarize_fam(path: Path) -> tuple[int, frozenset[tuple[str, str]], int, int, int, int, int, int]:
    sample_count = 0
    sample_ids: set[tuple[str, str]] = set()
    sex_1 = 0
    sex_2 = 0
    sex_other = 0
    pheno_1 = 0
    pheno_2 = 0
    pheno_other = 0

    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            sample_count += 1
            sample_ids.add((fields[0], fields[1]))

            sex = fields[4]
            if sex == "1":
                sex_1 += 1
            elif sex == "2":
                sex_2 += 1
            else:
                sex_other += 1

            pheno = fields[5]
            if pheno == "1":
                pheno_1 += 1
            elif pheno == "2":
                pheno_2 += 1
            else:
                pheno_other += 1

    return sample_count, frozenset(sample_ids), sex_1, sex_2, sex_other, pheno_1, pheno_2, pheno_other


def relative_display_path(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return os.path.relpath(path, base)


def collect_study_summary(analysis_root: Path, study_id: str) -> StudySummary:
    stage1_dir = analysis_root / study_id / "stage1"
    prefix = stage1_dir / study_id
    bed_path = prefix.with_suffix(".bed")
    bim_path = prefix.with_suffix(".bim")
    fam_path = prefix.with_suffix(".fam")
    summary_path = stage1_dir / "summary.txt"

    sample_count = None
    sample_ids: frozenset[tuple[str, str]] = frozenset()
    variant_count = None
    sex_1 = 0
    sex_2 = 0
    sex_other = 0
    pheno_1 = 0
    pheno_2 = 0
    pheno_other = 0

    if fam_path.exists():
        (
            sample_count,
            sample_ids,
            sex_1,
            sex_2,
            sex_other,
            pheno_1,
            pheno_2,
            pheno_other,
        ) = summarize_fam(fam_path)

    if bim_path.exists():
        variant_count = count_lines(bim_path)

    return StudySummary(
        study_id=study_id,
        bed_exists=bed_path.exists(),
        bim_exists=bim_path.exists(),
        fam_exists=fam_path.exists(),
        summary_exists=summary_path.exists(),
        sample_count=sample_count,
        sample_ids=sample_ids,
        variant_count=variant_count,
        sex_1=sex_1,
        sex_2=sex_2,
        sex_other=sex_other,
        pheno_1=pheno_1,
        pheno_2=pheno_2,
        pheno_other=pheno_other,
    )


def fmt_int(value: int | None) -> str:
    if value is None:
        return "-"
    return f"{value:,}"


def build_markdown(repo_root: Path, analysis_root: Path, summaries: list[StudySummary]) -> str:
    complete = [item for item in summaries if item.complete]
    incomplete = [item for item in summaries if not item.complete]

    total_samples = sum(item.sample_count or 0 for item in complete)
    unique_sample_ids: set[tuple[str, str]] = set()
    for item in summaries:
        unique_sample_ids.update(item.sample_ids)

    lines: list[str] = []
    lines.append("# Stage 1 Summary\n\n")
    lines.append(f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n")
    lines.append(f"Analysis root: `{relative_display_path(analysis_root, repo_root)}`\n\n")

    lines.append("## Overview\n\n")
    lines.append(f"- Expected studies: {len(summaries)}\n")
    lines.append(f"- Complete stage-1 outputs: {len(complete)}\n")
    lines.append(f"- Incomplete or missing studies: {len(incomplete)}\n")
    lines.append(f"- Total samples across complete studies: {total_samples:,}\n")
    lines.append(f"- Total unique samples across studies with final `.fam` outputs: {len(unique_sample_ids):,}\n\n")

    lines.append("## Complete Studies\n\n")
    if complete:
        lines.append(
            "| Study | Samples | Variants | Sex 1 | Sex 2 | Sex Other | Pheno 1 | Pheno 2 | Pheno Other |\n"
        )
        lines.append(
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n"
        )
        for item in sorted(complete, key=lambda value: value.study_id):
            lines.append(
                f"| {item.study_id} | {fmt_int(item.sample_count)} | {fmt_int(item.variant_count)} | "
                f"{item.sex_1:,} | {item.sex_2:,} | {item.sex_other:,} | "
                f"{item.pheno_1:,} | {item.pheno_2:,} | {item.pheno_other:,} |\n"
            )
        lines.append("\n")
    else:
        lines.append("No complete stage-1 outputs were found.\n\n")

    lines.append("## Incomplete Or Missing Studies\n\n")
    if incomplete:
        lines.append("| Study | Status | Samples | Variants |\n")
        lines.append("| --- | --- | ---: | ---: |\n")
        for item in sorted(incomplete, key=lambda value: value.study_id):
            lines.append(
                f"| {item.study_id} | {item.file_status} | {fmt_int(item.sample_count)} | {fmt_int(item.variant_count)} |\n"
            )
        lines.append("\n")
    else:
        lines.append("All expected stage-1 studies are complete.\n\n")

    lines.append("## Notes\n\n")
    lines.append("- `Samples` are counted from the final `.fam` file.\n")
    lines.append("- `Total unique samples` are counted from unique final `(FID, IID)` pairs across all available stage-1 `.fam` files.\n")
    lines.append("- `Variants` are counted from the final `.bim` file.\n")
    lines.append("- `Sex 1` and `Sex 2` come from column 5 of the final `.fam` file.\n")
    lines.append("- `Pheno 1` and `Pheno 2` come from column 6 of the final `.fam` file.\n")
    lines.append("- `Sex Other` and `Pheno Other` capture any value outside `1` or `2`.\n")
    return "".join(lines)


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    analysis_root = Path(args.analysis_root).resolve()
    scripts_dir = Path(args.scripts_dir).resolve()
    output_path = Path(args.output).resolve()

    study_ids = read_study_ids(scripts_dir)
    summaries = [collect_study_summary(analysis_root, study_id) for study_id in study_ids]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_markdown(repo_root, analysis_root, summaries))

    complete_count = sum(1 for item in summaries if item.complete)
    print(f"Wrote {relative_display_path(output_path, repo_root)}")
    print(f"Expected studies: {len(summaries)}")
    print(f"Complete studies: {complete_count}")
    print(f"Incomplete studies: {len(summaries) - complete_count}")


if __name__ == "__main__":
    main()
