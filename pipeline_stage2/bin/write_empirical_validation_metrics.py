#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
from pathlib import Path


VARIANT_FIELDNAMES = [
    "chrom",
    "pos",
    "variant_id",
    "maf",
    "n_samples",
    "empirical_dosage_r2",
    "dose0",
    "dose0_n",
    "mean_observed_dosage",
    "mean_imputed_dosage",
    "dosage_bias",
    "imputation_r2",
]
SUMMARY_FIELDNAMES = [
    "study",
    "chrom",
    "validation_variants",
    "empirical_r2_variants",
    "mean_empirical_dosage_r2",
    "median_empirical_dosage_r2",
    "p5_empirical_dosage_r2",
    "dose0_variants",
    "mean_dose0",
    "median_dose0",
    "mean_abs_dosage_bias",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare typed genotypes with imputed dosages at overlapping variants.")
    parser.add_argument("--target-vcf", required=True)
    parser.add_argument("--imputed-vcf", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--variant-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--min-samples", type=int, default=20)
    parser.add_argument("--min-dose0-samples", type=int, default=5)
    return parser.parse_args()


def open_vcf(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def fmt_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return ""
    return f"{value:.12g}"


def parse_info_float(info_field: str, keys: tuple[str, ...]) -> float | None:
    for entry in info_field.split(";"):
        if "=" not in entry:
            continue
        key, raw_value = entry.split("=", 1)
        if key not in keys:
            continue
        try:
            return float(raw_value)
        except ValueError:
            return None
    return None


def split_gt(gt_value: str) -> list[str]:
    if gt_value in {"", "."}:
        return []
    return gt_value.replace("|", "/").split("/")


def target_dosage(sample_field: str, format_keys: list[str], *, flipped: bool) -> tuple[float | None, int]:
    if "GT" not in format_keys:
        return None, 0
    values = sample_field.split(":")
    gt_index = format_keys.index("GT")
    if len(values) <= gt_index:
        return None, 0

    alleles = split_gt(values[gt_index])
    if not alleles or any(allele == "." for allele in alleles):
        return None, 0
    try:
        allele_counts = [int(allele) for allele in alleles]
    except ValueError:
        return None, 0
    if any(allele not in {0, 1} for allele in allele_counts):
        return None, 0

    ploidy = len(allele_counts)
    alt_count = sum(1 for allele in allele_counts if allele == 1)
    dosage = ploidy - alt_count if flipped else alt_count
    return float(dosage), ploidy


def imputed_dosage(sample_field: str, format_keys: list[str]) -> float | None:
    values = sample_field.split(":")
    for key in ("DS", "ADS"):
        if key in format_keys:
            index = format_keys.index(key)
            if len(values) > index and values[index] not in {"", "."}:
                try:
                    return float(values[index].split(",")[-1])
                except ValueError:
                    return None

    if "HDS" in format_keys:
        index = format_keys.index("HDS")
        if len(values) > index and values[index] not in {"", "."}:
            parts = values[index].replace("|", ",").split(",")
            try:
                return sum(float(part) for part in parts if part not in {"", "."})
            except ValueError:
                return None

    if "GP" in format_keys:
        index = format_keys.index("GP")
        if len(values) > index and values[index] not in {"", "."}:
            try:
                probs = [float(part) for part in values[index].replace("|", ",").split(",")]
            except ValueError:
                return None
            if len(probs) >= 3:
                return probs[1] + 2.0 * probs[2]
            if len(probs) == 2:
                return probs[1]

    return None


def correlation_squared(x_values: list[float], y_values: list[float]) -> float | None:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return None
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    ss_x = sum((value - x_mean) ** 2 for value in x_values)
    ss_y = sum((value - y_mean) ** 2 for value in y_values)
    if ss_x <= 0.0 or ss_y <= 0.0:
        return None
    cov = sum((x_value - x_mean) * (y_value - y_mean) for x_value, y_value in zip(x_values, y_values))
    r = cov / math.sqrt(ss_x * ss_y)
    return max(0.0, min(1.0, r * r))


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


class VcfIterator:
    def __init__(self, path: Path):
        self.path = path
        self.handle = open_vcf(path)
        self.samples: list[str] = []
        self.next_fields: list[str] | None = None
        self._read_header()

    def _read_header(self) -> None:
        for line in self.handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                self.samples = fields[9:]
                return
        raise RuntimeError(f"No #CHROM header found in {self.path}")

    def close(self) -> None:
        self.handle.close()

    def read_record(self) -> list[str] | None:
        if self.next_fields is not None:
            fields = self.next_fields
            self.next_fields = None
            return fields
        for line in self.handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9:
                return fields
        return None


def record_pos(fields: list[str]) -> int:
    return int(fields[1])


def load_target_group(iterator: VcfIterator, first_record: list[str] | None = None) -> tuple[int | None, list[list[str]]]:
    first = first_record if first_record is not None else iterator.read_record()
    if first is None:
        return None, []
    pos = record_pos(first)
    group = [first]
    while True:
        record = iterator.read_record()
        if record is None:
            return pos, group
        next_pos = record_pos(record)
        if next_pos != pos:
            iterator.next_fields = record
            return pos, group
        group.append(record)


def allele_match(target_fields: list[str], imputed_fields: list[str]) -> bool | None:
    target_ref = target_fields[3].upper()
    target_alt = target_fields[4].upper()
    imputed_ref = imputed_fields[3].upper()
    imputed_alt = imputed_fields[4].upper()
    if target_ref == imputed_ref and target_alt == imputed_alt:
        return False
    if target_ref == imputed_alt and target_alt == imputed_ref:
        return True
    return None


def write_validation_metrics(
    target_vcf: Path,
    imputed_vcf: Path,
    *,
    study: str,
    chrom: str,
    variant_output_path: Path,
    summary_output_path: Path,
    min_samples: int,
    min_dose0_samples: int,
) -> None:
    target = VcfIterator(target_vcf)
    imputed = VcfIterator(imputed_vcf)
    variant_rows: list[dict[str, object]] = []

    target_sample_index = {sample: index for index, sample in enumerate(target.samples)}
    sample_pairs = [(target_sample_index[sample], index) for index, sample in enumerate(imputed.samples) if sample in target_sample_index]
    if not sample_pairs:
        raise RuntimeError("Target and imputed VCFs have no overlapping sample IDs")

    target_group_pos, target_group = load_target_group(target)
    try:
        while True:
            imputed_fields = imputed.read_record()
            if imputed_fields is None:
                break
            imputed_pos = record_pos(imputed_fields)

            while target_group_pos is not None and target_group_pos < imputed_pos:
                target_group_pos, target_group = load_target_group(target)
            if target_group_pos != imputed_pos:
                continue

            match: tuple[list[str], bool] | None = None
            for target_fields in target_group:
                flipped = allele_match(target_fields, imputed_fields)
                if flipped is not None:
                    match = (target_fields, flipped)
                    break
            if match is None:
                continue

            target_fields, flipped = match
            target_format = target_fields[8].split(":")
            imputed_format = imputed_fields[8].split(":")
            observed: list[float] = []
            imputed_doses: list[float] = []
            total_alt_dosage = 0.0
            total_ploidy = 0

            for target_index, imputed_index in sample_pairs:
                target_sample = target_fields[9 + target_index]
                imputed_sample = imputed_fields[9 + imputed_index]
                observed_dose, ploidy = target_dosage(target_sample, target_format, flipped=flipped)
                predicted_dose = imputed_dosage(imputed_sample, imputed_format)
                if observed_dose is None or predicted_dose is None:
                    continue
                observed.append(observed_dose)
                imputed_doses.append(predicted_dose)
                total_alt_dosage += observed_dose
                total_ploidy += ploidy

            if not observed:
                continue

            allele_frequency = total_alt_dosage / total_ploidy if total_ploidy else None
            maf = min(allele_frequency, 1.0 - allele_frequency) if allele_frequency is not None else None
            r2 = correlation_squared(observed, imputed_doses) if len(observed) >= min_samples else None
            dose0_values = [predicted for truth, predicted in zip(observed, imputed_doses) if truth == 0.0]
            dose0 = sum(dose0_values) / len(dose0_values) if len(dose0_values) >= min_dose0_samples else None
            mean_observed = sum(observed) / len(observed)
            mean_imputed = sum(imputed_doses) / len(imputed_doses)
            variant_rows.append(
                {
                    "chrom": imputed_fields[0],
                    "pos": imputed_fields[1],
                    "variant_id": imputed_fields[2],
                    "maf": fmt_float(maf),
                    "n_samples": len(observed),
                    "empirical_dosage_r2": fmt_float(r2),
                    "dose0": fmt_float(dose0),
                    "dose0_n": len(dose0_values),
                    "mean_observed_dosage": fmt_float(mean_observed),
                    "mean_imputed_dosage": fmt_float(mean_imputed),
                    "dosage_bias": fmt_float(mean_imputed - mean_observed),
                    "imputation_r2": fmt_float(parse_info_float(imputed_fields[7], ("R2", "INFO", "DR2", "ER2"))),
                }
            )
    finally:
        target.close()
        imputed.close()

    variant_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    with variant_output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(variant_rows)

    empirical_values = [float(row["empirical_dosage_r2"]) for row in variant_rows if row["empirical_dosage_r2"] != ""]
    dose0_values = [float(row["dose0"]) for row in variant_rows if row["dose0"] != ""]
    abs_bias_values = [abs(float(row["dosage_bias"])) for row in variant_rows if row["dosage_bias"] != ""]
    summary_row = {
        "study": study,
        "chrom": chrom,
        "validation_variants": len(variant_rows),
        "empirical_r2_variants": len(empirical_values),
        "mean_empirical_dosage_r2": fmt_float(statistics.fmean(empirical_values) if empirical_values else None),
        "median_empirical_dosage_r2": fmt_float(statistics.median(empirical_values) if empirical_values else None),
        "p5_empirical_dosage_r2": fmt_float(percentile(empirical_values, 5)),
        "dose0_variants": len(dose0_values),
        "mean_dose0": fmt_float(statistics.fmean(dose0_values) if dose0_values else None),
        "median_dose0": fmt_float(statistics.median(dose0_values) if dose0_values else None),
        "mean_abs_dosage_bias": fmt_float(statistics.fmean(abs_bias_values) if abs_bias_values else None),
    }
    with summary_output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerow(summary_row)


def main() -> None:
    args = parse_args()
    write_validation_metrics(
        Path(args.target_vcf),
        Path(args.imputed_vcf),
        study=args.study,
        chrom=args.chrom,
        variant_output_path=Path(args.variant_output),
        summary_output_path=Path(args.summary_output),
        min_samples=args.min_samples,
        min_dose0_samples=args.min_dose0_samples,
    )


if __name__ == "__main__":
    main()
