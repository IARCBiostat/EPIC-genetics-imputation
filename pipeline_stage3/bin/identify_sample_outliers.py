#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Identify heterozygosity and ancestry outliers.")
    parser.add_argument("--eigenvec", required=True)
    parser.add_argument("--het", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--pc-count", type=int, default=10)
    parser.add_argument("--ancestry-z-threshold", type=float, default=6.0)
    parser.add_argument("--het-sd-threshold", type=float, default=3.0)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        header: list[str] | None = None
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if header is None:
                header = [value.lstrip("#") for value in parts]
                continue
            row = {header[idx]: parts[idx] if idx < len(parts) else "" for idx in range(len(header))}
            rows.append(row)
    return rows


def median(values: list[float]) -> float:
    sorted_values = sorted(values)
    n = len(sorted_values)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2 == 1:
        return sorted_values[mid]
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2.0


def stdev(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(variance, 0.0))


def robust_scale(values: list[float]) -> tuple[float, float]:
    center = median(values)
    deviations = [abs(value - center) for value in values]
    mad = median(deviations)
    if mad > 0:
        return center, 1.4826 * mad
    return center, stdev(values)


def write_id_file(path: Path, ids: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for fid, iid in ids:
            handle.write(f"{fid}\t{iid}\n")


def write_detail_file(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def identify_heterozygosity_outliers(
    rows: list[dict[str, str]],
    sd_threshold: float,
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    scored_rows: list[tuple[dict[str, str], float]] = []
    for row in rows:
        fid = row.get("FID", "")
        iid = row.get("IID", "")
        f_coeff = row.get("F", "")
        score: float | None = None

        if f_coeff not in {"", "NA"}:
            try:
                score = float(f_coeff)
            except ValueError:
                score = None

        if score is None:
            obs_key = "O(HOM)" if "O(HOM)" in row else "O.HOM."
            nm_key = "N(NM)" if "N(NM)" in row else "N.NM."
            try:
                obs_hom = float(row.get(obs_key, ""))
                n_nm = float(row.get(nm_key, ""))
                if n_nm > 0:
                    score = (n_nm - obs_hom) / n_nm
            except ValueError:
                score = None

        if score is not None:
            scored_rows.append(({"FID": fid, "IID": iid}, score))

    scores = [score for _, score in scored_rows]
    mean = sum(scores) / len(scores) if scores else 0.0
    sd = stdev(scores)

    outliers: list[tuple[str, str]] = []
    details: list[dict[str, str]] = []
    for sample, score in scored_rows:
        z_score = 0.0 if sd == 0 else (score - mean) / sd
        if sd > 0 and abs(z_score) > sd_threshold:
            outliers.append((sample["FID"], sample["IID"]))
            details.append(
                {
                    "FID": sample["FID"],
                    "IID": sample["IID"],
                    "HET_SCORE": f"{score:.6f}",
                    "Z_SCORE": f"{z_score:.6f}",
                }
            )

    return outliers, details


def identify_ancestry_outliers(
    rows: list[dict[str, str]],
    pc_count: int,
    z_threshold: float,
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    pc_names = [f"PC{index}" for index in range(1, pc_count + 1)]
    available_pcs = [pc for pc in pc_names if rows and pc in rows[0]]
    pc_values: dict[str, list[float]] = {pc: [] for pc in available_pcs}

    for row in rows:
        for pc in available_pcs:
            try:
                pc_values[pc].append(float(row[pc]))
            except ValueError:
                pc_values[pc].append(0.0)

    pc_centers_scales = {pc: robust_scale(values) for pc, values in pc_values.items()}

    outliers: list[tuple[str, str]] = []
    details: list[dict[str, str]] = []

    for row in rows:
        fid = row.get("FID", "")
        iid = row.get("IID", "")
        flagged_pcs: list[str] = []
        max_abs_z = 0.0

        for pc in available_pcs:
            try:
                value = float(row[pc])
            except ValueError:
                value = 0.0
            center, scale = pc_centers_scales[pc]
            z_score = 0.0 if scale == 0 else (value - center) / scale
            max_abs_z = max(max_abs_z, abs(z_score))
            if scale > 0 and abs(z_score) > z_threshold:
                flagged_pcs.append(pc)

        if flagged_pcs:
            outliers.append((fid, iid))
            details.append(
                {
                    "FID": fid,
                    "IID": iid,
                    "MAX_ABS_ROBUST_Z": f"{max_abs_z:.6f}",
                    "FLAGGED_PCS": ",".join(flagged_pcs),
                }
            )

    return outliers, details


def main() -> None:
    args = parse_args()
    out_prefix = Path(args.out_prefix)

    het_rows = read_table(Path(args.het)) if Path(args.het).exists() else []
    eigenvec_rows = read_table(Path(args.eigenvec)) if Path(args.eigenvec).exists() else []

    het_ids, het_details = identify_heterozygosity_outliers(het_rows, args.het_sd_threshold)
    ancestry_ids, ancestry_details = identify_ancestry_outliers(
        eigenvec_rows,
        args.pc_count,
        args.ancestry_z_threshold,
    )

    write_id_file(out_prefix.with_suffix(".heterozygosity_outliers.id"), het_ids)
    write_id_file(out_prefix.with_suffix(".ancestry_outliers.id"), ancestry_ids)

    write_detail_file(
        out_prefix.with_suffix(".heterozygosity_outliers.tsv"),
        ["FID", "IID", "HET_SCORE", "Z_SCORE"],
        het_details,
    )
    write_detail_file(
        out_prefix.with_suffix(".ancestry_outliers.tsv"),
        ["FID", "IID", "MAX_ABS_ROBUST_Z", "FLAGGED_PCS"],
        ancestry_details,
    )


if __name__ == "__main__":
    main()
