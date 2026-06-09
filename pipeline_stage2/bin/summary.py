#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parent
REPO_ROOT = BIN_DIR.parents[1]
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from common import parse_studies  # noqa: E402
EXPECTED_CHROMS = [str(chrom) for chrom in range(1, 23)] + ["X"]
QUALITY_KEYS = ("R2", "INFO", "DR2", "ER2")
RSID_RE = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
CHRPOS_RE = re.compile(r"^(?:chr)?[^:]+:\d+:[^:]+:[^:]+$")
ID_SAMPLE_LIMIT = 1000


@dataclass
class VcfScan:
    sample_count: int | None
    variant_count: int
    quality_tag_counts: dict[str, int]
    quality_scored_variants: int
    quality_sum: float
    quality_min: float | None
    quality_max: float | None
    quality_ge_min: int
    quality_ge_high: int
    id_style_counts: dict[str, int]


@dataclass
class StudySummary:
    study_id: str
    expected_chroms: list[str]
    stage1_samples: int | None
    stage2_samples: int | None
    chr_files: int
    chr_indexes: int
    total_variants: int
    total_size_bytes: int
    empty_chroms: list[str]
    missing_chroms: list[str]
    sample_mismatch: bool
    id_style: str
    quality_tag: str | None
    quality_scored_variants: int
    quality_ge_min: int
    quality_ge_high: int
    quality_mean: float | None
    quality_min: float | None
    quality_max: float | None

    @property
    def complete(self) -> bool:
        return (
            self.chr_files == len(self.expected_chroms)
            and self.chr_indexes == len(self.expected_chroms)
            and not self.empty_chroms
            and not self.missing_chroms
            and not self.sample_mismatch
            and self.stage2_samples is not None
            and self.quality_scored_variants > 0
        )

    @property
    def status(self) -> str:
        parts: list[str] = []
        if self.chr_files != len(self.expected_chroms):
            parts.append(f"vcf {self.chr_files}/{len(self.expected_chroms)}")
        if self.chr_indexes != len(self.expected_chroms):
            parts.append(f"index {self.chr_indexes}/{len(self.expected_chroms)}")
        if self.empty_chroms:
            parts.append(f"empty: {','.join(self.empty_chroms)}")
        if self.missing_chroms:
            parts.append(f"missing: {','.join(self.missing_chroms)}")
        if self.sample_mismatch:
            parts.append("sample-mismatch")
        if self.quality_scored_variants == 0 and self.total_variants > 0:
            parts.append("quality-missing")
        if not parts:
            return "complete"
        return "; ".join(parts)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write a markdown summary of stage-2 imputation outputs.")
    parser.add_argument(
        "--analysis-root",
        default=str(REPO_ROOT / "analysis"),
        help="Root directory containing analysis/<STUDY>/stage1 and stage2 outputs.",
    )
    parser.add_argument(
        "--stage1-root",
        default=None,
        help="Optional root directory containing stage-1 outputs. Defaults to --analysis-root.",
    )
    parser.add_argument(
        "--output",
        default=str(REPO_ROOT / "analysis" / "stage2-summary.md"),
        help="Path to the markdown summary to write.",
    )
    parser.add_argument(
        "--studies",
        default="all",
        help="Comma-separated study IDs to summarize. Default: all studies found under the stage-1 root.",
    )
    parser.add_argument(
        "--chromosomes",
        default="all",
        help="Comma-separated chromosomes to summarize. Default: all chromosomes present in stage 1.",
    )
    parser.add_argument(
        "--min-r2-threshold",
        type=float,
        default=0.3,
        help="Primary imputation quality threshold for pass/fail counts.",
    )
    parser.add_argument(
        "--high-quality-threshold",
        type=float,
        default=0.8,
        help="High-confidence imputation quality threshold.",
    )
    return parser.parse_args()


def project_relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def fmt_int(value: int | None) -> str:
    if value is None:
        return "-"
    return f"{value:,}"


def fmt_size_gib(value: int) -> str:
    return f"{value / (1024 ** 3):.2f}"


def fmt_float(value: float | None, digits: int = 3) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def fmt_percent(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "-"
    return f"{(100.0 * numerator / denominator):.1f}%"


def parse_study_ids(analysis_root: Path, study_arg: str) -> list[str]:
    return parse_studies(study_arg, analysis_root)


def parse_chromosomes(chromosome_arg: str) -> list[str] | None:
    if chromosome_arg == "all":
        return None
    chroms = []
    for item in chromosome_arg.split(","):
        chrom = item.strip().replace("chr", "")
        if not chrom:
            continue
        chrom = "X" if chrom == "23" else chrom
        if chrom not in EXPECTED_CHROMS:
            raise SystemExit(f"Unsupported chromosome for summary: {item}")
        chroms.append(chrom)
    return chroms


def count_stage1_samples(fam_path: Path) -> int | None:
    if not fam_path.exists():
        return None
    with fam_path.open() as handle:
        return sum(1 for line in handle if line.strip())


def detect_stage1_chroms(bim_path: Path) -> list[str]:
    if not bim_path.exists():
        return EXPECTED_CHROMS.copy()

    observed: set[str] = set()
    with bim_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom = line.split(maxsplit=1)[0]
            if chrom in {str(value) for value in range(1, 23)}:
                observed.add(chrom)
            elif chrom in {"X", "23", "chrX"}:
                observed.add("X")

    expected = [chrom for chrom in EXPECTED_CHROMS if chrom in observed]
    return expected or EXPECTED_CHROMS.copy()


def count_vcf_samples(vcf_path: Path) -> int | None:
    if not vcf_path.exists():
        return None
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return max(0, len(fields) - 9)
    return None


def classify_id(value: str) -> str:
    if value in {"", "."}:
        return "missing"
    if RSID_RE.match(value):
        return "rsID"
    if CHRPOS_RE.match(value):
        return "chr:pos:REF:ALT"
    return "other"


def parse_quality(info_field: str) -> tuple[str | None, float | None]:
    values: dict[str, float | None] = {}
    for entry in info_field.split(";"):
        if "=" not in entry:
            continue
        key, raw_value = entry.split("=", 1)
        if key not in QUALITY_KEYS:
            continue
        try:
            values[key] = float(raw_value)
        except ValueError:
            values[key] = None

    for key in QUALITY_KEYS:
        if key in values:
            return key, values[key]
    return None, None


def detect_id_style(id_style_counts: dict[str, int]) -> str:
    chr_style = id_style_counts.get("chr:pos:REF:ALT", 0)
    rs_style = id_style_counts.get("rsID", 0)
    other_style = id_style_counts.get("other", 0)
    missing_style = id_style_counts.get("missing", 0)

    informative = chr_style + rs_style + other_style
    if informative == 0 and missing_style > 0:
        return "missing"
    if chr_style > 0 and rs_style == 0 and other_style == 0:
        return "chr:pos:REF:ALT"
    if rs_style > 0 and chr_style == 0 and other_style == 0:
        return "rsID"
    if informative == 0:
        return "-"
    return "mixed"


def dominant_quality_tag(quality_tag_counts: dict[str, int]) -> str | None:
    non_zero = {key: value for key, value in quality_tag_counts.items() if value > 0}
    if not non_zero:
        return None
    return sorted(non_zero.items(), key=lambda item: (-item[1], item[0]))[0][0]


def scan_vcf(
    vcf_path: Path,
    min_r2_threshold: float,
    high_quality_threshold: float,
    id_sample_limit: int,
) -> VcfScan:
    sample_count: int | None = None
    variant_count = 0
    quality_tag_counts = {key: 0 for key in QUALITY_KEYS}
    quality_scored_variants = 0
    quality_sum = 0.0
    quality_min: float | None = None
    quality_max: float | None = None
    quality_ge_min = 0
    quality_ge_high = 0
    id_style_counts = {"rsID": 0, "chr:pos:REF:ALT": 0, "other": 0, "missing": 0}
    sampled_ids = 0

    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample_count = max(0, len(fields) - 9)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) < 8:
                continue

            variant_count += 1

            if sampled_ids < id_sample_limit:
                id_style_counts[classify_id(fields[2])] += 1
                sampled_ids += 1

            quality_tag, quality_value = parse_quality(fields[7])
            if quality_tag:
                quality_tag_counts[quality_tag] += 1
            if quality_value is None:
                continue

            quality_scored_variants += 1
            quality_sum += quality_value
            if quality_min is None or quality_value < quality_min:
                quality_min = quality_value
            if quality_max is None or quality_value > quality_max:
                quality_max = quality_value
            if quality_value >= min_r2_threshold:
                quality_ge_min += 1
            if quality_value >= high_quality_threshold:
                quality_ge_high += 1

    return VcfScan(
        sample_count=sample_count,
        variant_count=variant_count,
        quality_tag_counts=quality_tag_counts,
        quality_scored_variants=quality_scored_variants,
        quality_sum=quality_sum,
        quality_min=quality_min,
        quality_max=quality_max,
        quality_ge_min=quality_ge_min,
        quality_ge_high=quality_ge_high,
        id_style_counts=id_style_counts,
    )


def int_from_metrics(value: str | None, default: int = 0) -> int:
    if value in {None, ""}:
        return default
    return int(value)


def optional_int_from_metrics(value: str | None) -> int | None:
    if value in {None, ""}:
        return None
    return int(value)


def float_from_metrics(value: str | None) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def metrics_threshold_matches(value: str | None, expected: float) -> bool:
    if value in {None, ""}:
        return False
    try:
        return abs(float(value) - expected) < 1e-12
    except ValueError:
        return False


def read_staged_metrics(
    metrics_path: Path,
    min_r2_threshold: float,
    high_quality_threshold: float,
) -> VcfScan | None:
    if not metrics_path.exists():
        return None

    with metrics_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)

    if not row:
        return None

    if not metrics_threshold_matches(row.get("min_r2_threshold"), min_r2_threshold):
        return None
    if not metrics_threshold_matches(row.get("high_quality_threshold"), high_quality_threshold):
        return None

    try:
        return VcfScan(
            sample_count=optional_int_from_metrics(row.get("sample_count")),
            variant_count=int_from_metrics(row.get("variant_count")),
            quality_tag_counts={
                "R2": int_from_metrics(row.get("quality_tag_R2")),
                "INFO": int_from_metrics(row.get("quality_tag_INFO")),
                "DR2": int_from_metrics(row.get("quality_tag_DR2")),
                "ER2": int_from_metrics(row.get("quality_tag_ER2")),
            },
            quality_scored_variants=int_from_metrics(row.get("quality_scored_variants")),
            quality_sum=float(row.get("quality_sum") or 0.0),
            quality_min=float_from_metrics(row.get("quality_min")),
            quality_max=float_from_metrics(row.get("quality_max")),
            quality_ge_min=int_from_metrics(row.get("quality_ge_min")),
            quality_ge_high=int_from_metrics(row.get("quality_ge_high")),
            id_style_counts={
                "rsID": int_from_metrics(row.get("id_rsID")),
                "chr:pos:REF:ALT": int_from_metrics(row.get("id_chrpos")),
                "other": int_from_metrics(row.get("id_other")),
                "missing": int_from_metrics(row.get("id_missing")),
            },
        )
    except ValueError:
        return None


def collect_study_summary(
    analysis_root: Path,
    stage1_root: Path,
    study_id: str,
    min_r2_threshold: float,
    high_quality_threshold: float,
    selected_chroms: list[str] | None,
) -> StudySummary:
    stage1_fam = stage1_root / study_id / "stage1" / f"{study_id}.fam"
    stage1_bim = stage1_root / study_id / "stage1" / f"{study_id}.bim"
    stage2_dir = analysis_root / study_id / "stage2"
    expected_chroms = detect_stage1_chroms(stage1_bim)
    if selected_chroms is not None:
        expected_chroms = [chrom for chrom in expected_chroms if chrom in selected_chroms]

    stage1_samples = count_stage1_samples(stage1_fam)
    stage2_samples: int | None = None
    chr_files = 0
    chr_indexes = 0
    total_variants = 0
    total_size_bytes = 0
    empty_chroms: list[str] = []
    missing_chroms: list[str] = []
    quality_tag_counts = {key: 0 for key in QUALITY_KEYS}
    quality_scored_variants = 0
    quality_sum = 0.0
    quality_min: float | None = None
    quality_max: float | None = None
    quality_ge_min = 0
    quality_ge_high = 0
    id_style_counts = {"rsID": 0, "chr:pos:REF:ALT": 0, "other": 0, "missing": 0}
    remaining_id_samples = ID_SAMPLE_LIMIT

    for chrom in expected_chroms:
        vcf_path = stage2_dir / f"{study_id}_chr{chrom}_GxS.imputed.vcf.gz"
        tbi_path = Path(f"{vcf_path}.tbi")

        if not vcf_path.exists():
            missing_chroms.append(chrom)
            continue

        chr_files += 1
        total_size_bytes += vcf_path.stat().st_size

        if tbi_path.exists():
            chr_indexes += 1

        metrics_path = stage2_dir / "report" / "tables" / f"chr{chrom}.imputation_metrics.tsv"
        legacy_metrics_path = stage2_dir / "report" / "tables" / f"{study_id}_chr{chrom}.imputation_metrics.tsv"
        if not metrics_path.exists() and legacy_metrics_path.exists():
            metrics_path = legacy_metrics_path
        scan = read_staged_metrics(
            metrics_path,
            min_r2_threshold=min_r2_threshold,
            high_quality_threshold=high_quality_threshold,
        )
        if scan is None:
            scan = scan_vcf(
                vcf_path,
                min_r2_threshold=min_r2_threshold,
                high_quality_threshold=high_quality_threshold,
                id_sample_limit=remaining_id_samples,
            )

        if stage2_samples is None:
            stage2_samples = scan.sample_count

        total_variants += scan.variant_count
        if scan.variant_count == 0:
            empty_chroms.append(chrom)

        for key, value in scan.quality_tag_counts.items():
            quality_tag_counts[key] += value
        quality_scored_variants += scan.quality_scored_variants
        quality_sum += scan.quality_sum
        if scan.quality_min is not None and (quality_min is None or scan.quality_min < quality_min):
            quality_min = scan.quality_min
        if scan.quality_max is not None and (quality_max is None or scan.quality_max > quality_max):
            quality_max = scan.quality_max
        quality_ge_min += scan.quality_ge_min
        quality_ge_high += scan.quality_ge_high

        for key, value in scan.id_style_counts.items():
            id_style_counts[key] += value
        remaining_id_samples = max(0, remaining_id_samples - sum(scan.id_style_counts.values()))

    sample_mismatch = (
        stage1_samples is not None
        and stage2_samples is not None
        and stage1_samples != stage2_samples
    )
    quality_mean = quality_sum / quality_scored_variants if quality_scored_variants else None

    return StudySummary(
        study_id=study_id,
        expected_chroms=expected_chroms,
        stage1_samples=stage1_samples,
        stage2_samples=stage2_samples,
        chr_files=chr_files,
        chr_indexes=chr_indexes,
        total_variants=total_variants,
        total_size_bytes=total_size_bytes,
        empty_chroms=empty_chroms,
        missing_chroms=missing_chroms,
        sample_mismatch=sample_mismatch,
        id_style=detect_id_style(id_style_counts),
        quality_tag=dominant_quality_tag(quality_tag_counts),
        quality_scored_variants=quality_scored_variants,
        quality_ge_min=quality_ge_min,
        quality_ge_high=quality_ge_high,
        quality_mean=quality_mean,
        quality_min=quality_min,
        quality_max=quality_max,
    )


def build_markdown(
    analysis_root: Path,
    stage1_root: Path,
    summaries: list[StudySummary],
    min_r2_threshold: float,
    high_quality_threshold: float,
) -> str:
    complete = [item for item in summaries if item.complete]
    incomplete = [item for item in summaries if not item.complete]

    total_variants = sum(item.total_variants for item in complete)
    total_size_bytes = sum(item.total_size_bytes for item in complete)
    total_quality_scored = sum(item.quality_scored_variants for item in complete)
    total_quality_ge_min = sum(item.quality_ge_min for item in complete)
    total_quality_ge_high = sum(item.quality_ge_high for item in complete)

    lines: list[str] = []
    lines.append("# Stage 2 Summary\n\n")
    lines.append(f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}\n\n")
    lines.append(f"Analysis root: `{project_relative(analysis_root)}`\n\n")
    lines.append(f"Stage-1 root: `{project_relative(stage1_root)}`\n\n")

    lines.append("## Overview\n\n")
    lines.append(f"- Expected studies: {len(summaries):,}\n")
    lines.append(f"- Complete stage-2 outputs: {len(complete):,}\n")
    lines.append(f"- Incomplete or missing studies: {len(incomplete):,}\n")
    lines.append(f"- Total imputed variants across complete studies: {total_variants:,}\n")
    lines.append(
        f"- Variants meeting quality threshold (`>= {min_r2_threshold}`) across complete studies: "
        f"{total_quality_ge_min:,} / {total_quality_scored:,} ({fmt_percent(total_quality_ge_min, total_quality_scored)})\n"
    )
    lines.append(
        f"- Variants meeting high-quality threshold (`>= {high_quality_threshold}`) across complete studies: "
        f"{total_quality_ge_high:,} / {total_quality_scored:,} ({fmt_percent(total_quality_ge_high, total_quality_scored)})\n"
    )
    lines.append(f"- Total stage-2 size across complete studies: {fmt_size_gib(total_size_bytes)} GiB\n\n")

    lines.append("## Complete Studies\n\n")
    if complete:
        lines.append(
            "| Study | Stage1 Samples | Stage2 Samples | Sample Match | ID Style | Qual Tag | Variants | >= "
            f"{min_r2_threshold} | >= {high_quality_threshold} | Mean Qual | Chr Files | Indexed | Size GiB |\n"
        )
        lines.append("| --- | ---: | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for item in sorted(complete, key=lambda value: value.study_id):
            lines.append(
                f"| {item.study_id} | {fmt_int(item.stage1_samples)} | {fmt_int(item.stage2_samples)} | "
                f"{'yes' if not item.sample_mismatch else 'no'} | "
                f"{item.id_style} | {item.quality_tag or '-'} | "
                f"{item.total_variants:,} | "
                f"{item.quality_ge_min:,} ({fmt_percent(item.quality_ge_min, item.quality_scored_variants)}) | "
                f"{item.quality_ge_high:,} ({fmt_percent(item.quality_ge_high, item.quality_scored_variants)}) | "
                f"{fmt_float(item.quality_mean)} | "
                f"{item.chr_files}/{len(item.expected_chroms)} | {item.chr_indexes}/{len(item.expected_chroms)} | "
                f"{fmt_size_gib(item.total_size_bytes)} |\n"
            )
        lines.append("\n")
    else:
        lines.append("No complete stage-2 outputs were found.\n\n")

    lines.append("## Incomplete Or Missing Studies\n\n")
    if incomplete:
        lines.append(
            "| Study | Status | Stage1 Samples | Stage2 Samples | ID Style | Qual Tag | Variants | >= "
            f"{min_r2_threshold} | Mean Qual | Chr Files | Indexed | Empty Chrs | Missing Chrs |\n"
        )
        lines.append("| --- | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |\n")
        for item in sorted(incomplete, key=lambda value: value.study_id):
            lines.append(
                f"| {item.study_id} | {item.status} | {fmt_int(item.stage1_samples)} | {fmt_int(item.stage2_samples)} | "
                f"{item.id_style} | {item.quality_tag or '-'} | "
                f"{item.total_variants:,} | {item.quality_ge_min:,} ({fmt_percent(item.quality_ge_min, item.quality_scored_variants)}) | "
                f"{fmt_float(item.quality_mean)} | "
                f"{item.chr_files}/{len(item.expected_chroms)} | {item.chr_indexes}/{len(item.expected_chroms)} | "
                f"{','.join(item.empty_chroms) if item.empty_chroms else '-'} | "
                f"{','.join(item.missing_chroms) if item.missing_chroms else '-'} |\n"
            )
        lines.append("\n")
    else:
        lines.append("All expected stage-2 studies are complete.\n\n")

    lines.append("## Notes\n\n")
    lines.append("- `Stage1 Samples` are counted from `analysis/<STUDY>/stage1/<STUDY>.fam`.\n")
    lines.append("- `Stage2 Samples` are counted from the sample columns in the imputed VCF header.\n")
    lines.append("- `ID Style` is inferred from sampled final VCF IDs and reports `rsID`, `chr:pos:REF:ALT`, or `mixed`.\n")
    lines.append(
        "- `Qual Tag` is the INFO field used for imputation quality, prioritising `R2`, then `INFO`, `DR2`, and `ER2`.\n"
    )
    lines.append("- `Chr Files` expects the chromosomes present in the stage-1 `.bim`; studies with no chrX in stage 1 are assessed on autosomes only.\n")
    lines.append("- `Indexed` counts `.tbi` files present beside the final imputed VCFs.\n")
    lines.append("- `Variants` are summed across the final imputed stage-2 VCFs.\n")
    lines.append("- Per-chromosome counts are read from `analysis/<STUDY>/stage2/report/tables/*imputation_metrics.tsv` when present; missing or threshold-mismatched metrics fall back to a direct VCF scan.\n")
    lines.append(f"- `>= {min_r2_threshold}` counts variants meeting the pipeline-quality threshold.\n")
    lines.append(f"- `>= {high_quality_threshold}` counts variants in a higher-confidence quality tier.\n")
    lines.append("- `Mean Qual` is the mean imputation-quality value across scored variants.\n")
    lines.append("- `Empty Chrs` flags any output VCF with zero variant records.\n")
    lines.append("- `Missing Chrs` flags chromosomes with no final imputed VCF.\n")
    lines.append("- `Sample Match` should stay `yes`; imputation should not change the sample count.\n")
    return "".join(lines)


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    stage1_root = Path(args.stage1_root).resolve() if args.stage1_root else analysis_root
    output_path = Path(args.output).resolve()

    study_ids = parse_study_ids(analysis_root, args.studies)
    selected_chroms = parse_chromosomes(args.chromosomes)
    summaries = [
        collect_study_summary(
            analysis_root,
            stage1_root,
            study_id,
            min_r2_threshold=args.min_r2_threshold,
            high_quality_threshold=args.high_quality_threshold,
            selected_chroms=selected_chroms,
        )
        for study_id in study_ids
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        build_markdown(
            analysis_root,
            stage1_root,
            summaries,
            min_r2_threshold=args.min_r2_threshold,
            high_quality_threshold=args.high_quality_threshold,
        )
    )

    complete_count = sum(1 for item in summaries if item.complete)
    print(f"Wrote {output_path}")
    print(f"Expected studies: {len(summaries)}")
    print(f"Complete studies: {complete_count}")
    print(f"Incomplete studies: {len(summaries) - complete_count}")


if __name__ == "__main__":
    main()
