#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_CHROMS = [str(chrom) for chrom in range(1, 23)] + ["X"]
RSID_RE = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
CHRPOS_RE = re.compile(r"^(?:chr)?[^:]+:\d+:[^:]+:[^:]+$")
STAGE2_VCF_RE = re.compile(r"^.+_chr([0-9X]+)_GxS\.imputed\.vcf\.gz$")


@dataclass
class StudySummary:
    study_id: str
    stage1_samples: int | None
    stage3_samples: int | None
    final_variants: int | None
    input_variants: int
    post_r2_maf_variants: int
    hwe_exclude_variants: int
    rsid_variants: int
    fallback_variants: int
    duplicate_rsid_fallbacks: int
    related_identified: int
    related_removed: int
    heterozygosity_outliers: int
    ancestry_identified: int
    ancestry_removed: int
    total_removed: int
    variant_qc_files: int
    hwe_no_controls: bool
    complete: bool
    status: str


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Write a markdown summary of stage-3 outputs.")
    parser.add_argument("--analysis-root", default=str(repo_root / "analysis"))
    parser.add_argument("--stage1-root", default=None)
    parser.add_argument("--stage2-root", default=None)
    parser.add_argument("--studies", default="all")
    parser.add_argument("--output", default=str(repo_root / "analysis" / "stage3-summary.md"))
    return parser.parse_args()


def fmt_int(value: int | None) -> str:
    if value is None:
        return "-"
    return f"{value:,}"


def count_non_header_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def count_fam_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def get_stage2_chroms(stage2_root: Path, study_id: str) -> list[str]:
    """Return sorted chromosomes present in stage2 VCFs for this study, or EXPECTED_CHROMS if none found."""
    stage2_dir = stage2_root / study_id / "stage2"
    if not stage2_dir.exists():
        return []
    chroms: set[str] = set()
    for vcf in stage2_dir.iterdir():
        m = STAGE2_VCF_RE.match(vcf.name)
        if m:
            chroms.add(m.group(1))
    if not chroms:
        return []

    def _sort_key(c: str) -> int:
        return 23 if c == "X" else int(c)

    return sorted(chroms, key=_sort_key)


def parse_study_ids(analysis_root: Path, stage1_root: Path, stage2_root: Path, studies_arg: str) -> list[str]:
    if studies_arg != "all":
        return sorted({item.strip() for item in studies_arg.split(",") if item.strip()})

    # Stage 1 is the canonical entry point; only studies processed there are expected.
    # Require the per-study .fam file to exist — this excludes cohort reporting dirs
    # (analysis/cohort/stage1/) which are created by stage1 reporting but are not studies.
    return sorted(
        path.parent.name
        for path in stage1_root.glob("*/stage1")
        if path.is_dir() and (path / f"{path.parent.name}.fam").exists()
    )


def load_tsv_row(path: Path) -> dict[str, str] | None:
    if not path.exists():
        return None
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return row
    return None


def collect_study_summary(analysis_root: Path, stage1_root: Path, stage2_root: Path, study_id: str) -> StudySummary:
    stage1_fam = stage1_root / study_id / "stage1" / f"{study_id}.fam"
    stage3_dir = analysis_root / study_id / "stage3"
    final_dir = stage3_dir / "final"

    study_expected_chroms = get_stage2_chroms(stage2_root, study_id) or EXPECTED_CHROMS

    report_tables_dir = stage3_dir / "report" / "tables"
    variant_qc_files = sorted(report_tables_dir.glob("*.variant_metrics.tsv"))
    sample_qc_path = report_tables_dir / f"{study_id}.sample_review.tsv"
    sample_qc_row = load_tsv_row(sample_qc_path) or {}

    input_variants = 0
    post_r2_maf_variants = 0
    hwe_exclude_variants = 0
    duplicate_rsid_fallbacks = 0
    rsid_variants = 0
    fallback_variants = 0
    hwe_no_controls = False

    for path in variant_qc_files:
        row = load_tsv_row(path)
        if not row:
            continue
        input_variants += int(row["input_variants"])
        post_r2_maf_variants += int(row["post_r2_maf_variants"])
        hwe_exclude_variants += int(row.get("hwe_exclude_count", 0))
        duplicate_rsid_fallbacks += int(row["duplicate_rsid_fallbacks"])
        rsid_variants += int(row["rsid_variants"])
        fallback_variants += int(row["fallback_variants"])
        if row.get("hwe_no_controls", "0") == "1":
            hwe_no_controls = True

    stage1_samples = count_fam_lines(stage1_fam)

    pre_final_samples = int(sample_qc_row.get("pre_final_samples", 0) or 0)
    total_removed_count = int(sample_qc_row.get("total_removed", 0) or 0)
    stage3_samples = max(pre_final_samples - total_removed_count, 0) if pre_final_samples else None

    final_variants = post_r2_maf_variants if post_r2_maf_variants else None

    related_identified = int(sample_qc_row.get("related_identified", 0) or 0)
    related_removed = int(sample_qc_row.get("related_removed", 0) or 0)
    heterozygosity_outliers = int(sample_qc_row.get("heterozygosity_outliers", 0) or 0)
    ancestry_identified = int(sample_qc_row.get("ancestry_outliers_identified", 0) or 0)
    ancestry_removed = int(sample_qc_row.get("ancestry_outliers_removed", 0) or 0)
    total_removed = int(sample_qc_row.get("total_removed", 0) or 0)

    final_pgen_count = len(list(final_dir.glob("*.pgen"))) if final_dir.exists() else 0

    status_parts: list[str] = []
    if final_pgen_count != len(study_expected_chroms):
        status_parts.append(f"final-pgen {final_pgen_count}/{len(study_expected_chroms)}")
    if len(variant_qc_files) != len(study_expected_chroms):
        status_parts.append(f"variant-metrics {len(variant_qc_files)}/{len(study_expected_chroms)}")
    if not sample_qc_path.exists():
        status_parts.append("sample-review-missing")
    if stage1_samples is not None and stage3_samples is not None and stage3_samples > stage1_samples:
        status_parts.append("sample-count-increase")

    complete = not status_parts
    status = "complete" if complete else "; ".join(status_parts)

    return StudySummary(
        study_id=study_id,
        stage1_samples=stage1_samples,
        stage3_samples=stage3_samples,
        final_variants=final_variants,
        input_variants=input_variants,
        post_r2_maf_variants=post_r2_maf_variants,
        hwe_exclude_variants=hwe_exclude_variants,
        rsid_variants=rsid_variants,
        fallback_variants=fallback_variants,
        duplicate_rsid_fallbacks=duplicate_rsid_fallbacks,
        related_identified=related_identified,
        related_removed=related_removed,
        heterozygosity_outliers=heterozygosity_outliers,
        ancestry_identified=ancestry_identified,
        ancestry_removed=ancestry_removed,
        total_removed=total_removed,
        variant_qc_files=len(variant_qc_files),
        hwe_no_controls=hwe_no_controls,
        complete=complete,
        status=status,
    )


def build_markdown(analysis_root: Path, summaries: list[StudySummary]) -> str:
    complete = [item for item in summaries if item.complete]
    incomplete = [item for item in summaries if not item.complete]

    total_final_samples = sum(item.stage3_samples or 0 for item in complete)
    total_final_variants = sum(item.final_variants or 0 for item in complete)

    lines: list[str] = []
    lines.append("# Stage 3 Summary\n\n")
    lines.append(f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n")
    lines.append(f"Analysis root: `{analysis_root}`\n\n")

    lines.append("## Overview\n\n")
    lines.append(f"- Expected studies: {len(summaries)}\n")
    lines.append(f"- Complete stage-3 outputs: {len(complete)}\n")
    lines.append(f"- Incomplete or missing studies: {len(incomplete)}\n")
    lines.append(f"- Total final samples across complete studies: {total_final_samples:,}\n")
    lines.append(f"- Total final variants across complete studies: {total_final_variants:,}\n")
    hwe_skipped = [item for item in summaries if item.hwe_no_controls]
    if hwe_skipped:
        skipped_names = ", ".join(item.study_id for item in sorted(hwe_skipped, key=lambda x: x.study_id))
        lines.append(f"- **Note — HWE exclusion list empty (no controls):** {skipped_names}\n")
    lines.append("\n")

    lines.append("## Complete Studies\n\n")
    if complete:
        lines.append(
            "| Study | Stage1 Samples | Stage3 Samples | Variants In | Post R2/MAF | HWE Excl | Final Variants | rsID | Fallback ID | Related ID | Related Removed | Het | Ancestry ID | Ancestry Removed | Total Removed | HWE Skip |\n"
        )
        lines.append(
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |\n"
        )
        for item in sorted(complete, key=lambda value: value.study_id):
            hwe_skip_flag = "no controls" if item.hwe_no_controls else ""
            lines.append(
                f"| {item.study_id} | {fmt_int(item.stage1_samples)} | {fmt_int(item.stage3_samples)} | "
                f"{fmt_int(item.input_variants)} | {fmt_int(item.post_r2_maf_variants)} | {fmt_int(item.hwe_exclude_variants)} | "
                f"{fmt_int(item.final_variants)} | {fmt_int(item.rsid_variants)} | {fmt_int(item.fallback_variants)} | "
                f"{item.related_identified:,} | {item.related_removed:,} | {item.heterozygosity_outliers:,} | "
                f"{item.ancestry_identified:,} | {item.ancestry_removed:,} | {item.total_removed:,} | {hwe_skip_flag} |\n"
            )
        lines.append("\n")
    else:
        lines.append("No complete stage-3 outputs were found.\n\n")

    lines.append("## Incomplete Or Missing Studies\n\n")
    if incomplete:
        lines.append("| Study | Status | Stage1 Samples | Stage3 Samples | Final Variants |\n")
        lines.append("| --- | --- | ---: | ---: | ---: |\n")
        for item in sorted(incomplete, key=lambda value: value.study_id):
            lines.append(
                f"| {item.study_id} | {item.status} | {fmt_int(item.stage1_samples)} | {fmt_int(item.stage3_samples)} | {fmt_int(item.final_variants)} |\n"
            )
        lines.append("\n")
    else:
        lines.append("All expected stage-3 studies are complete.\n\n")

    lines.append("## Notes\n\n")
    lines.append("- `Variants In` is the total count of stage-2 variants entering stage-3 QC.\n")
    lines.append("- `Post R2/MAF` is after imputation-quality and MAF filtering.\n")
    lines.append("- `HWE Excl` is the count of variants written to the per-study `hwe.exclude` file (flagged by autosomal HWE test in controls only); no hard filter is applied — variants pass through and the list is for downstream use.\n")
    lines.append("- `HWE Skip` is flagged `no controls` for studies where HWE was requested but the stage-1 FAM contains no control samples (phenotype=1); the exclusion list will be empty for these studies.\n")
    lines.append("- `Final Variants` is the sum of post-R2/MAF variants across all chromosomes (HWE does not reduce this count).\n")
    lines.append("- `rsID` counts final stage-3 variants with a retained dbSNP rsID.\n")
    lines.append("- `Fallback ID` counts final stage-3 variants using `chr:pos:ref:alt` identifiers.\n")
    lines.append("- `Ancestry ID` is always reported; `Ancestry Removed` is only non-zero when ancestry exclusion is enabled.\n")
    lines.append("- Completeness is determined by the presence of all per-chromosome `.pgen` files in `stage3/final/` and all chromosome variant-metrics files.\n")
    return "".join(lines)


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    stage1_root = Path(args.stage1_root).resolve() if args.stage1_root else analysis_root
    stage2_root = Path(args.stage2_root).resolve() if args.stage2_root else analysis_root
    output_path = Path(args.output).resolve()

    study_ids = parse_study_ids(analysis_root, stage1_root, stage2_root, args.studies)
    summaries = [collect_study_summary(analysis_root, stage1_root, stage2_root, study_id) for study_id in study_ids]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_markdown(analysis_root, summaries))

    complete_count = sum(1 for item in summaries if item.complete)
    print(f"Wrote {output_path}")
    print(f"Expected studies: {len(summaries)}")
    print(f"Complete studies: {complete_count}")
    print(f"Incomplete studies: {len(summaries) - complete_count}")


if __name__ == "__main__":
    main()
