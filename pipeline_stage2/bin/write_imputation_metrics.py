#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


QUALITY_KEYS = ("R2", "INFO", "DR2", "ER2")
RSID_RE = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
CHRPOS_RE = re.compile(r"^(?:chr)?[^:]+:\d+:[^:]+:[^:]+$")
ID_SAMPLE_LIMIT = 1000

FIELDNAMES = [
    "study",
    "chrom",
    "vcf_name",
    "index_present",
    "size_bytes",
    "sample_count",
    "variant_count",
    "quality_tag_R2",
    "quality_tag_INFO",
    "quality_tag_DR2",
    "quality_tag_ER2",
    "quality_scored_variants",
    "quality_sum",
    "quality_min",
    "quality_max",
    "quality_ge_min",
    "quality_ge_high",
    "min_r2_threshold",
    "high_quality_threshold",
    "id_rsID",
    "id_chrpos",
    "id_other",
    "id_missing",
]

VARIANT_FIELDNAMES = ["chrom", "pos", "variant_id", "r2", "maf", "empirical_dosage_r2", "dose0"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write staged metrics for one imputed VCF.")
    parser.add_argument("--vcf", required=True, help="Final imputed VCF, optionally bgzipped.")
    parser.add_argument("--study", required=True, help="Study identifier.")
    parser.add_argument("--chrom", required=True, help="Chromosome label.")
    parser.add_argument("--output", required=True, help="Output metrics TSV.")
    parser.add_argument("--variant-output", default=None, help="Optional per-variant R2/MAF TSV.")
    parser.add_argument("--min-r2-threshold", type=float, default=0.3)
    parser.add_argument("--high-quality-threshold", type=float, default=0.8)
    parser.add_argument("--id-sample-limit", type=int, default=ID_SAMPLE_LIMIT)
    return parser.parse_args()


def classify_id(value: str) -> str:
    if value in {"", "."}:
        return "missing"
    if RSID_RE.match(value):
        return "rsID"
    if CHRPOS_RE.match(value):
        return "chrpos"
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


def parse_info_float(info_field: str, target_key: str) -> float | None:
    for entry in info_field.split(";"):
        if "=" not in entry:
            continue
        key, raw_value = entry.split("=", 1)
        if key != target_key:
            continue
        try:
            return float(raw_value)
        except ValueError:
            return None
    return None


def parse_first_info_float(info_field: str, target_keys: tuple[str, ...]) -> float | None:
    for target_key in target_keys:
        value = parse_info_float(info_field, target_key)
        if value is not None:
            return value
    return None


def open_vcf(vcf_path: Path):
    if vcf_path.name.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    return vcf_path.open()


def scan_vcf_metrics(
    vcf_path: Path,
    *,
    study: str,
    chrom: str,
    min_r2_threshold: float,
    high_quality_threshold: float,
    id_sample_limit: int = ID_SAMPLE_LIMIT,
    variant_output_path: Path | None = None,
) -> dict[str, str | int | float]:
    sample_count: int | None = None
    variant_count = 0
    quality_tag_counts = {key: 0 for key in QUALITY_KEYS}
    quality_scored_variants = 0
    quality_sum = 0.0
    quality_min: float | None = None
    quality_max: float | None = None
    quality_ge_min = 0
    quality_ge_high = 0
    id_style_counts = {"rsID": 0, "chrpos": 0, "other": 0, "missing": 0}
    sampled_ids = 0

    variant_handle = None
    variant_writer = None
    if variant_output_path is not None:
        variant_output_path.parent.mkdir(parents=True, exist_ok=True)
        variant_handle = variant_output_path.open("w", newline="")
        variant_writer = csv.DictWriter(variant_handle, fieldnames=VARIANT_FIELDNAMES, delimiter="\t")
        variant_writer.writeheader()

    try:
        with open_vcf(vcf_path) as handle:
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
                maf_value = parse_info_float(fields[7], "MAF")
                empirical_value = parse_first_info_float(fields[7], ("ER2", "EMPIRICAL_R2", "EMP_R2", "LOO_R2"))
                dose0_value = parse_first_info_float(fields[7], ("Dose0", "DOSE0", "dose0", "DOSE_0"))
                if variant_writer is not None:
                    variant_writer.writerow(
                        {
                            "chrom": fields[0],
                            "pos": fields[1],
                            "variant_id": fields[2],
                            "r2": "" if quality_value is None else f"{quality_value:.12g}",
                            "maf": "" if maf_value is None else f"{maf_value:.12g}",
                            "empirical_dosage_r2": "" if empirical_value is None else f"{empirical_value:.12g}",
                            "dose0": "" if dose0_value is None else f"{dose0_value:.12g}",
                        }
                    )
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
    finally:
        if variant_handle is not None:
            variant_handle.close()

    return {
        "study": study,
        "chrom": chrom,
        "vcf_name": vcf_path.name,
        "index_present": int(Path(f"{vcf_path}.tbi").exists()),
        "size_bytes": vcf_path.stat().st_size,
        "sample_count": "" if sample_count is None else sample_count,
        "variant_count": variant_count,
        "quality_tag_R2": quality_tag_counts["R2"],
        "quality_tag_INFO": quality_tag_counts["INFO"],
        "quality_tag_DR2": quality_tag_counts["DR2"],
        "quality_tag_ER2": quality_tag_counts["ER2"],
        "quality_scored_variants": quality_scored_variants,
        "quality_sum": f"{quality_sum:.12g}",
        "quality_min": "" if quality_min is None else f"{quality_min:.12g}",
        "quality_max": "" if quality_max is None else f"{quality_max:.12g}",
        "quality_ge_min": quality_ge_min,
        "quality_ge_high": quality_ge_high,
        "min_r2_threshold": f"{min_r2_threshold:.12g}",
        "high_quality_threshold": f"{high_quality_threshold:.12g}",
        "id_rsID": id_style_counts["rsID"],
        "id_chrpos": id_style_counts["chrpos"],
        "id_other": id_style_counts["other"],
        "id_missing": id_style_counts["missing"],
    }


def write_metrics(
    vcf_path: Path,
    *,
    study: str,
    chrom: str,
    output_path: Path,
    min_r2_threshold: float,
    high_quality_threshold: float,
    id_sample_limit: int = ID_SAMPLE_LIMIT,
    variant_output_path: Path | None = None,
) -> None:
    row = scan_vcf_metrics(
        vcf_path,
        study=study,
        chrom=chrom,
        min_r2_threshold=min_r2_threshold,
        high_quality_threshold=high_quality_threshold,
        id_sample_limit=id_sample_limit,
        variant_output_path=variant_output_path,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    args = parse_args()
    write_metrics(
        Path(args.vcf),
        study=args.study,
        chrom=args.chrom,
        output_path=Path(args.output),
        min_r2_threshold=args.min_r2_threshold,
        high_quality_threshold=args.high_quality_threshold,
        id_sample_limit=args.id_sample_limit,
        variant_output_path=Path(args.variant_output) if args.variant_output else None,
    )


if __name__ == "__main__":
    main()
