#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


PHASING_INFO_KEYS = ("PQ", "GP", "PP", "PS")
VARIANT_FIELDNAMES = ["chrom", "pos", "variant_id", "metric_type", "metric_value"]
SUMMARY_FIELDNAMES = [
    "study",
    "chrom",
    "metric_type",
    "variant_count",
    "scored_variants",
    "mean_conf",
    "p5_conf",
    "n_low_conf_variants",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write staged phasing metrics for one phased VCF.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--variant-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def open_vcf(vcf_path: Path):
    if vcf_path.name.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    return vcf_path.open()


def parse_info_value(info_field: str, target_key: str) -> float | None:
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


def percentile(values: list[float], percent: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    rank = (len(ordered) - 1) * percent / 100.0
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - rank) + ordered[hi] * (rank - lo)


def fmt_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def write_phasing_metrics(
    vcf_path: Path,
    *,
    study: str,
    chrom: str,
    variant_output_path: Path,
    summary_output_path: Path,
) -> None:
    info_metric_key: str | None = None
    format_keys: list[str] = []
    first_sample_index: int | None = None
    previous_ps: str | None = None
    metric_type = "missing"
    variant_count = 0
    metric_values: list[float] = []

    variant_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)

    with variant_output_path.open("w", newline="") as variant_handle:
        variant_writer = csv.DictWriter(variant_handle, fieldnames=VARIANT_FIELDNAMES, delimiter="\t")
        variant_writer.writeheader()

        with open_vcf(vcf_path) as handle:
            for line in handle:
                if line.startswith("##INFO="):
                    for key in PHASING_INFO_KEYS:
                        if f"ID={key}," in line:
                            info_metric_key = info_metric_key or key
                    continue
                if line.startswith("#CHROM"):
                    header_fields = line.rstrip("\n").split("\t")
                    first_sample_index = 9 if len(header_fields) > 9 else None
                    continue
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue

                variant_count += 1
                metric_value: float | None = None
                if info_metric_key is not None:
                    metric_type = info_metric_key
                    metric_value = parse_info_value(fields[7], info_metric_key)
                elif first_sample_index is not None and len(fields) > first_sample_index:
                    format_keys = fields[8].split(":")
                    if "PS" in format_keys:
                        metric_type = "ps_switch_proxy"
                        ps_index = format_keys.index("PS")
                        sample_values = fields[first_sample_index].split(":")
                        ps_value = sample_values[ps_index] if len(sample_values) > ps_index else "."
                        metric_value = 0.0 if previous_ps is None or ps_value == previous_ps else 1.0
                        previous_ps = ps_value

                if metric_value is not None:
                    metric_values.append(metric_value)
                variant_writer.writerow(
                    {
                        "chrom": fields[0],
                        "pos": fields[1],
                        "variant_id": fields[2],
                        "metric_type": metric_type,
                        "metric_value": fmt_float(metric_value),
                    }
                )

    mean_conf = sum(metric_values) / len(metric_values) if metric_values else None
    p5_conf = percentile(metric_values, 5)
    n_low_conf = sum(1 for value in metric_values if value < 0.5)
    with summary_output_path.open("w", newline="") as summary_handle:
        writer = csv.DictWriter(summary_handle, fieldnames=SUMMARY_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "study": study,
                "chrom": chrom,
                "metric_type": metric_type,
                "variant_count": variant_count,
                "scored_variants": len(metric_values),
                "mean_conf": fmt_float(mean_conf),
                "p5_conf": fmt_float(p5_conf),
                "n_low_conf_variants": n_low_conf,
            }
        )


def main() -> None:
    args = parse_args()
    write_phasing_metrics(
        Path(args.vcf),
        study=args.study,
        chrom=args.chrom,
        variant_output_path=Path(args.variant_output),
        summary_output_path=Path(args.summary_output),
    )


if __name__ == "__main__":
    main()
