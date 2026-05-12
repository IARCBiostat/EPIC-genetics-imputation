#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import math
import re
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = Path(__file__).resolve().parent
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from common import (  # noqa: E402
    ASSET_COLUMNS,
    FLAG_COLUMNS,
    SPECIAL_STUDY_NOTES,
    asset_row,
    check_dependencies,
    ensure_dir,
    find_reference_vcf,
    flag_row,
    flags_output_path,
    manifest_output_path,
    merge_stage_study_assets,
    merge_stage_study_flags,
    parse_studies,
    read_reference_af_rows,
    run_command,
    safe_float,
    save_figure,
    setup_logger,
    stage_dir,
    study_qc_dirs,
    write_tsv,
)


COLUMN_RENAME_MAP = {
    "maf_bin": "MAF",
    "maf_midpoint": "MAF midpoint",
    "n_variants": "N variants",
    "scored_variants": "N variants",
    "n_total_variants": "Total variants",
    "n_r2_ge_0.3": "N variants R2 >= 0.3",
    "n_r2_ge_0.8": "N variants R2 >= 0.8",
    "n_empirical_r2": "N variants empirical R2",
    "n_dose0": "N variants Dose0",
    "dose0_variants": "N variants Dose0",
    "empirical_r2_variants": "N variants empirical R2",
    "validation_variants": "N variants",
    "mean_empirical_dosage_r2": "Mean empirical R2",
    "median_empirical_dosage_r2": "Median empirical R2",
    "mean_dose0": "Mean Dose0",
    "median_dose0": "Median Dose0",
    "mean_abs_dosage_bias": "Mean |dosage bias|",
    "mean_r2": "Mean R2",
    "median_r2": "Median R2",
    "p5_r2": "P5 R2",
    "mean_conf": "Mean confidence",
    "p5_conf": "P5 confidence",
    "n_low_conf_variants": "N variants low confidence",
    "imputation_basis": "N variants basis",
    "imputation_target": "N variants target",
    "total": "N variants",
    "chromosome": "Chromosome",
    "metric_type": "Metric",
    "imputed_af": "Imputed AF",
    "kg_af": "1000G AF",
    "af_diff": "AF diff",
    "pearson_r": "Pearson r",
    "slope": "Slope",
    "n_outlier_variants": "N outlier variants",
    "empirical_dosage_r2_mean_median": "Empirical R2 (mean/median)",
    "dose0_mean_median": "Dose0 (mean/median)",
    "imputation_r2_mean_median": "Imputation R2 (mean/median)",
}

plt.rcParams["figure.dpi"] = 150
sns.set_theme(style="whitegrid")

STAGE_NAME = "stage2"
PIPELINE_NAME = "pipeline_stage2"
TASK_SEQUENCE = [
    "variant_summary",
    "r2_by_maf",
    "af_concordance",
]
STAGE2_PY_PACKAGES = ["numpy", "pandas", "matplotlib", "seaborn"]
CHROMS = [str(chrom) for chrom in range(1, 23)] + ["X"]
R2_BINS = [0.0, 0.3, 0.5, 0.8, 0.9, 1.0]
R2_BIN_LABELS = ["bin_0_0.3", "bin_0.3_0.5", "bin_0.5_0.8", "bin_0.8_0.9", "bin_0.9_1.0"]
MAF_BIN_WIDTH = 0.01
MAF_BINS = [round(value * MAF_BIN_WIDTH, 2) for value in range(0, 51)]


def format_numeric(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, (int, np.integer)):
        return f"{value:,}"
    try:
        # Check if it's an integer-like float
        fval = float(value)
        if fval.is_integer():
            return f"{int(fval):,}"
        if fval == 0:
            return "0.0000"
        if abs(fval) < 0.0001:
            return f"{fval:.4E}"
        return f"{fval:.4f}"
    except (TypeError, ValueError):
        return str(value)


def _rename_headers(headers: list[str]) -> list[str]:
    return [COLUMN_RENAME_MAP.get(h, h) for h in headers]
ACTIVE_CHROMS: set[str] | None = None
STALE_LABEL_RE = re.compile(r"(?:^|_)[CD][0-9](?:_|\.|$)")
OBSOLETE_REPORT_FILES = {
    "genotyped_r2_by_maf.png",
    "r2_miami.png",
    "r2_distribution.png",
    "r2_windows.png",
    "validation_summary.tsv",
    "phasing_quality.assets.tsv",
    "phasing_plots.assets.tsv",
}
OBSOLETE_REPORT_TASKS = {
    "genotyped_r2_by_maf",
    "r2_distribution",
    "r2_windows",
    "validation_summary",
    "phasing_quality",
    "phasing_plots",
}


def _parse_chromosomes(chromosomes_arg: str) -> set[str] | None:
    if chromosomes_arg == "all":
        return None
    chroms: set[str] = set()
    for item in chromosomes_arg.split(","):
        chrom = item.strip().replace("chr", "")
        if chrom:
            chroms.add("X" if chrom == "23" else chrom)
    return chroms


def _chrom_from_path(path: Path) -> str | None:
    matcher = re.search(r"(?:^|_)chr([0-9X]+)", path.name)
    return matcher.group(1) if matcher else None


def _chrom_allowed(path: Path) -> bool:
    chrom = _chrom_from_path(path)
    return ACTIVE_CHROMS is None or chrom in ACTIVE_CHROMS


def _normalise_chrom(value: object) -> str:
    chrom = str(value).replace("chr", "").strip()
    return "X" if chrom == "23" else chrom


def _chrom_sort_key(chrom: object) -> tuple[int, str]:
    chrom_clean = _normalise_chrom(chrom)
    if chrom_clean == "X":
        return 23, chrom_clean
    try:
        return int(chrom_clean), chrom_clean
    except ValueError:
        return 99, chrom_clean


def _study_ids(studies_arg: str, analysis_root: Path) -> list[str]:
    return parse_studies(studies_arg, analysis_root)


def _project_relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _report_dirs(analysis_root: Path, study: str) -> dict[str, Path]:
    dirs = study_qc_dirs(analysis_root, study, STAGE_NAME)
    ensure_dir(dirs["tables"] / "imputation")
    ensure_dir(dirs["tables"] / "phasing")
    ensure_dir(dirs["tables"] / "imputation_validation")
    ensure_dir(dirs["figures"] / "imputation")
    ensure_dir(dirs["figures"] / "phasing")
    return dirs


def _replace_or_remove(source: Path, target: Path) -> None:
    if not source.exists() or source == target:
        return
    ensure_dir(target.parent)
    if target.exists():
        source.unlink()
        return
    source.replace(target)


def _normalise_staged_report_inputs(analysis_root: Path, study: str) -> None:
    report = stage_dir(analysis_root, study, STAGE_NAME) / "report"
    if not report.exists():
        return

    tables = report / "tables"
    for path in sorted(tables.glob(f"{study}_chr*.imputation_metrics.tsv")):
        chrom = _chrom_from_path(path)
        if chrom:
            _replace_or_remove(path, tables / f"chr{chrom}.imputation_metrics.tsv")

    for path in sorted((tables / "imputation").glob(f"{study}_chr*.r2_maf.tsv")):
        chrom = _chrom_from_path(path)
        if chrom:
            _replace_or_remove(path, tables / "imputation" / f"chr{chrom}.r2_maf.tsv")

    for path in sorted((tables / "phasing").glob(f"{study}_chr*_phase_quality.tsv")):
        chrom = _chrom_from_path(path)
        if chrom:
            _replace_or_remove(path, tables / "phasing" / f"chr{chrom}.phase_quality.tsv")

    for path in sorted((tables / "phasing").glob(f"{study}_chr*_phasing_metrics.tsv")):
        chrom = _chrom_from_path(path)
        if chrom:
            _replace_or_remove(path, tables / "phasing" / f"chr{chrom}.phasing_metrics.tsv")


def _cleanup_stale_report_outputs(analysis_root: Path, study: str) -> None:
    report = stage_dir(analysis_root, study, STAGE_NAME) / "report"
    if not report.exists():
        return
    for path in sorted(report.rglob("*"), reverse=True):
        if not path.is_file():
            continue
        task_name = path.name.rsplit(".", 2)[0]
        if (
            path.name == ".DS_Store"
            or path.suffix.lower() == ".pdf"
            or study in path.name
            or STALE_LABEL_RE.search(path.name)
            or path.name in OBSOLETE_REPORT_FILES
            or task_name in OBSOLETE_REPORT_TASKS
        ):
            path.unlink()

    _migrate_legacy_reference_cache(analysis_root)
    cohort_report = analysis_root / "cohort" / STAGE_NAME / "report"
    if cohort_report.exists():
        shutil.rmtree(cohort_report)


def _migrate_legacy_reference_cache(analysis_root: Path) -> None:
    legacy_cache_dir = analysis_root / "cohort" / STAGE_NAME / "report" / "tables" / "reference_cache"
    if not legacy_cache_dir.exists():
        return

    cache_dir = ensure_dir(analysis_root / "cohort" / STAGE_NAME / "reference" / "af_cache")
    for legacy_cache_path in sorted(legacy_cache_dir.glob("1kg_chr*_af.tsv.gz")):
        cache_path = cache_dir / legacy_cache_path.name
        if not cache_path.exists():
            shutil.copy2(legacy_cache_path, cache_path)


def _prepare_study_report_dir(analysis_root: Path, study: str) -> dict[str, Path]:
    _normalise_staged_report_inputs(analysis_root, study)
    _cleanup_stale_report_outputs(analysis_root, study)
    return _report_dirs(analysis_root, study)


def _note_flags(study: str, task: str) -> list[dict[str, object]]:
    if study not in SPECIAL_STUDY_NOTES:
        return []
    return [
        flag_row(
            study=study,
            task=task,
            metric="study_note",
            value=study,
            threshold="special-case metadata",
            flag_level="WARN",
            message=SPECIAL_STUDY_NOTES[study],
        )
    ]


def _stage2_imputation_variant_tables(analysis_root: Path, study: str) -> list[Path]:
    table_dir = stage_dir(analysis_root, study, STAGE_NAME) / "report" / "tables" / "imputation"
    return [path for path in sorted(table_dir.glob("chr*.r2_maf.tsv"), key=lambda item: _chrom_sort_key(_chrom_from_path(item) or "")) if _chrom_allowed(path)]


def _stage2_phasing_metric_tables(analysis_root: Path, study: str) -> list[Path]:
    table_dir = stage_dir(analysis_root, study, STAGE_NAME) / "report" / "tables" / "phasing"
    return [path for path in sorted(table_dir.glob("chr*.phasing_metrics.tsv"), key=lambda item: _chrom_sort_key(_chrom_from_path(item) or "")) if _chrom_allowed(path)]


def _stage2_phasing_variant_tables(analysis_root: Path, study: str) -> list[Path]:
    table_dir = stage_dir(analysis_root, study, STAGE_NAME) / "report" / "tables" / "phasing"
    return [path for path in sorted(table_dir.glob("chr*.phase_quality.tsv"), key=lambda item: _chrom_sort_key(_chrom_from_path(item) or "")) if _chrom_allowed(path)]


def _stage2_validation_variant_tables(analysis_root: Path, study: str) -> list[Path]:
    table_dir = stage_dir(analysis_root, study, STAGE_NAME) / "report" / "tables" / "imputation_validation"
    return [path for path in sorted(table_dir.glob("chr*.empirical_metrics.tsv"), key=lambda item: _chrom_sort_key(_chrom_from_path(item) or "")) if _chrom_allowed(path)]


def _load_or_build_r2_table(analysis_root: Path, study: str) -> Path | None:
    dirs = _report_dirs(analysis_root, study)
    output_path = dirs["tables"] / "r2_by_variant.tsv"

    staged_tables = _stage2_imputation_variant_tables(analysis_root, study)
    frames = [pd.read_csv(path, sep="\t") for path in staged_tables if path.stat().st_size]
    if frames:
        frame = pd.concat(frames, ignore_index=True)
        frame.to_csv(output_path, sep="\t", index=False)
        return output_path

    if output_path.exists() and ACTIVE_CHROMS is None:
        return output_path

    return None


def _load_or_build_empirical_table(analysis_root: Path, study: str) -> Path | None:
    dirs = _report_dirs(analysis_root, study)
    output_path = dirs["tables"] / "empirical_validation_by_variant.tsv"

    staged_tables = _stage2_validation_variant_tables(analysis_root, study)
    frames = [pd.read_csv(path, sep="\t") for path in staged_tables if path.stat().st_size]
    if frames:
        frame = pd.concat(frames, ignore_index=True)
        frame.to_csv(output_path, sep="\t", index=False)
        return output_path

    if output_path.exists() and ACTIVE_CHROMS is None:
        return output_path

    return None


def _write_manifest(dirs: dict[str, Path], task: str, assets: list[dict[str, object]]) -> None:
    write_tsv(assets, manifest_output_path(dirs, task), ASSET_COLUMNS)


def _write_flags(dirs: dict[str, Path], task: str, flags: list[dict[str, object]]) -> None:
    write_tsv(flags, flags_output_path(dirs, task), FLAG_COLUMNS)


def _asset(
    *,
    study: str,
    task: str,
    section: str,
    asset_id: str,
    asset_type: str,
    fmt: str,
    path: Path,
    title: str,
    caption: str,
    sort_order: int,
) -> dict[str, object]:
    row = asset_row(
        study=study,
        stage=STAGE_NAME,
        task=task,
        section=section,
        asset_id=asset_id,
        asset_type=asset_type,
        fmt=fmt,
        path=path,
        title=title,
        caption=caption,
        sort_order=sort_order,
    )
    row["path"] = _project_relative(path)
    return row


def _stage1_variant_sets(analysis_root: Path, study: str) -> tuple[set[str], set[str]]:
    bim_path = analysis_root / study / "stage1" / f"{study}.bim"
    allele_keys: set[str] = set()
    pos_keys: set[str] = set()
    if not bim_path.exists():
        return allele_keys, pos_keys

    with bim_path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            chrom = _normalise_chrom(fields[0])
            pos = fields[3]
            a1 = fields[4].upper()
            a2 = fields[5].upper()
            if len(a1) == 1 and len(a2) == 1:
                alleles = sorted([a1, a2])
                allele_keys.add(f"{chrom}:{pos}:{alleles[0]}:{alleles[1]}")
            pos_keys.add(f"{chrom}:{pos}")
    return allele_keys, pos_keys


def _variant_keys(row: pd.Series) -> tuple[str, str]:
    chrom = _normalise_chrom(row.get("chrom", ""))
    pos_val = row.get("pos_numeric", row.get("pos", ""))
    try:
        pos = str(int(pd.to_numeric(pos_val)))
    except (ValueError, TypeError):
        pos = str(pos_val)

    # Robustly handle CHROM:POS:REF:ALT or variants thereof
    variant_id = str(row.get("variant_id", ""))
    parts = variant_id.replace("chr", "").split(":")
    if len(parts) >= 4:
        chrom = _normalise_chrom(parts[0])
        pos = parts[1]
        alleles = sorted([parts[-2].upper(), parts[-1].upper()])
        return f"{chrom}:{pos}:{alleles[0]}:{alleles[1]}", f"{chrom}:{pos}"

    # Fallback to REF/ALT columns if available
    ref = str(row.get("ref", row.get("REF", ""))).upper()
    alt = str(row.get("alt", row.get("ALT", ""))).upper()
    if ref and alt:
        alleles = sorted([ref, alt])
        return f"{chrom}:{pos}:{alleles[0]}:{alleles[1]}", f"{chrom}:{pos}"

    return "", f"{chrom}:{pos}"


def _variant_merge_key(row: pd.Series) -> str:
    allele_key, pos_key = _variant_keys(row)
    return allele_key or pos_key


def _load_variant_frame(analysis_root: Path, study: str) -> tuple[pd.DataFrame, Path | None]:
    table_path = _load_or_build_r2_table(analysis_root, study)
    if table_path is None:
        return pd.DataFrame(), None

    frame = pd.read_csv(table_path, sep="\t")
    if frame.empty:
        return frame, table_path

    frame["r2"] = pd.to_numeric(frame["r2"] if "r2" in frame else pd.Series(dtype=float), errors="coerce")
    frame["maf"] = pd.to_numeric(frame["maf"] if "maf" in frame else pd.Series(dtype=float), errors="coerce")
    frame["pos_numeric"] = pd.to_numeric(frame["pos"] if "pos" in frame else pd.Series(dtype=float), errors="coerce")
    chrom_values = frame["chrom"] if "chrom" in frame else pd.Series([""] * len(frame), index=frame.index)
    frame["chrom_clean"] = chrom_values.astype(str).map(_normalise_chrom)
    allele_keys, pos_keys = _stage1_variant_sets(analysis_root, study)
    if allele_keys or pos_keys:
        keys = frame.apply(_variant_keys, axis=1)
        frame["is_genotyped"] = [allele_key in allele_keys or pos_key in pos_keys for allele_key, pos_key in keys]
    else:
        frame["is_genotyped"] = False

    validation_path = _load_or_build_empirical_table(analysis_root, study)
    if validation_path is not None:
        validation = pd.read_csv(validation_path, sep="\t")
        if not validation.empty:
            validation["pos_numeric"] = pd.to_numeric(validation["pos"] if "pos" in validation else pd.Series(dtype=float), errors="coerce")
            validation["validation_key"] = validation.apply(_variant_merge_key, axis=1)
            validation_columns = [
                "validation_key",
                "empirical_dosage_r2",
                "dose0",
                "dose0_n",
                "n_samples",
                "mean_observed_dosage",
                "mean_imputed_dosage",
                "dosage_bias",
                "imputation_r2",
            ]
            available_columns = [column for column in validation_columns if column in validation]
            for column in available_columns:
                if column != "validation_key":
                    validation[column] = pd.to_numeric(validation[column], errors="coerce")

            frame["validation_key"] = frame.apply(_variant_merge_key, axis=1)
            
            # Drop existing validation columns from frame to avoid name collisions (_x, _y)
            cols_to_drop = [c for c in available_columns if c in frame and c != "validation_key"]
            if cols_to_drop:
                frame.drop(columns=cols_to_drop, inplace=True)

            frame = frame.merge(
                validation[available_columns].drop_duplicates("validation_key"),
                on="validation_key",
                how="left",
            )
            validation_hit = pd.Series(False, index=frame.index)
            for column in ("empirical_dosage_r2", "dose0", "n_samples"):
                if column in frame:
                    validation_hit = validation_hit | frame[column].notna()
            frame["is_genotyped"] = frame["is_genotyped"] | validation_hit
            frame.drop(columns=["validation_key"], inplace=True)

    frame["is_imputed"] = ~frame["is_genotyped"]
    if table_path is not None and ACTIVE_CHROMS is None:
        frame.to_csv(table_path, sep="\t", index=False)
    return frame, table_path


def _phasing_basis_counts(analysis_root: Path, study: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for table in _stage2_phasing_metric_tables(analysis_root, study):
        if not table.exists() or table.stat().st_size == 0:
            continue
        frame = pd.read_csv(table, sep="\t")
        for record in frame.to_dict("records"):
            chrom = _normalise_chrom(record.get("chrom", _chrom_from_path(table) or ""))
            try:
                counts[chrom] = counts.get(chrom, 0) + int(record.get("variant_count", 0) or 0)
            except (TypeError, ValueError):
                continue
    return counts


def _maf_bin_summary(frame: pd.DataFrame, value_column: str = "r2") -> pd.DataFrame:
    data = frame.dropna(subset=["maf", value_column]).copy()
    if data.empty:
        return pd.DataFrame(columns=["maf_bin", "maf_midpoint", "n_variants", f"mean_{value_column}"])
    data["maf_clipped"] = data["maf"].clip(0, 0.499999)
    data["maf_bin"] = pd.cut(data["maf_clipped"], bins=MAF_BINS, include_lowest=True, right=False)
    summary = (
        data.groupby("maf_bin", observed=True)
        .agg(n_variants=("variant_id", "size"), mean_value=(value_column, "mean"))
        .reset_index()
    )
    summary["maf_midpoint"] = summary["maf_bin"].map(lambda interval: float((interval.left + interval.right) / 2))
    summary["maf_bin"] = summary["maf_bin"].map(lambda interval: f"{interval.left:.2f}-{interval.right:.2f}")
    return summary.rename(columns={"mean_value": f"mean_{value_column}"})


def _empirical_maf_bin_summary(frame: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "maf_bin",
        "maf_midpoint",
        "n_variants",
        "n_empirical_r2",
        "mean_empirical_dosage_r2",
        "median_empirical_dosage_r2",
        "n_dose0",
        "mean_dose0",
        "median_dose0",
    ]
    data = frame.dropna(subset=["maf"]).copy()
    if data.empty:
        return pd.DataFrame(columns=columns)
    data["maf_clipped"] = data["maf"].clip(0, 0.499999)
    data["maf_bin_interval"] = pd.cut(data["maf_clipped"], bins=MAF_BINS, include_lowest=True, right=False)
    rows: list[dict[str, object]] = []
    for interval, subset in data.groupby("maf_bin_interval", observed=True):
        empirical = pd.to_numeric(subset["empirical_dosage_r2"], errors="coerce").dropna() if "empirical_dosage_r2" in subset else pd.Series(dtype=float)
        dose0 = pd.to_numeric(subset["dose0"], errors="coerce").dropna() if "dose0" in subset else pd.Series(dtype=float)
        rows.append(
            {
                "maf_bin": f"{interval.left:.2f}-{interval.right:.2f}",
                "maf_midpoint": float((interval.left + interval.right) / 2),
                "n_variants": len(subset),
                "n_empirical_r2": len(empirical),
                "mean_empirical_dosage_r2": float(empirical.mean()) if not empirical.empty else float("nan"),
                "median_empirical_dosage_r2": float(empirical.median()) if not empirical.empty else float("nan"),
                "n_dose0": len(dose0),
                "mean_dose0": float(dose0.mean()) if not dose0.empty else float("nan"),
                "median_dose0": float(dose0.median()) if not dose0.empty else float("nan"),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _empirical_validation_summary_rows(frame: pd.DataFrame) -> list[dict[str, object]]:
    columns = [
        "chromosome",
        "empirical_r2_variants",
        "mean_empirical_dosage_r2",
        "median_empirical_dosage_r2",
        "dose0_variants",
        "mean_dose0",
        "median_dose0",
        "mean_abs_dosage_bias",
    ]
    if frame.empty:
        return []
    data = frame.copy()
    for column in ("empirical_dosage_r2", "dose0", "dosage_bias"):
        if column in data:
            data[column] = pd.to_numeric(data[column], errors="coerce")
    rows: list[dict[str, object]] = []
    for chrom, subset in sorted(data.groupby("chrom_clean"), key=lambda item: _chrom_sort_key(item[0])):
        empirical = subset["empirical_dosage_r2"].dropna() if "empirical_dosage_r2" in subset else pd.Series(dtype=float)
        dose0 = subset["dose0"].dropna() if "dose0" in subset else pd.Series(dtype=float)
        bias = subset["dosage_bias"].dropna().abs() if "dosage_bias" in subset else pd.Series(dtype=float)
        rows.append(
            {
                "chromosome": chrom,
                "validation_variants": len(subset),
                "empirical_r2_variants": len(empirical),
                "mean_empirical_dosage_r2": float(empirical.mean()) if not empirical.empty else float("nan"),
                "median_empirical_dosage_r2": float(empirical.median()) if not empirical.empty else float("nan"),
                "dose0_variants": len(dose0),
                "mean_dose0": float(dose0.mean()) if not dose0.empty else float("nan"),
                "median_dose0": float(dose0.median()) if not dose0.empty else float("nan"),
                "mean_abs_dosage_bias": float(bias.mean()) if not bias.empty else float("nan"),
            }
        )
    if rows:
        empirical = data["empirical_dosage_r2"].dropna() if "empirical_dosage_r2" in data else pd.Series(dtype=float)
        dose0 = data["dose0"].dropna() if "dose0" in data else pd.Series(dtype=float)
        bias = data["dosage_bias"].dropna().abs() if "dosage_bias" in data else pd.Series(dtype=float)
        rows.append(
            {
                "chromosome": "Total",
                "validation_variants": len(data),
                "empirical_r2_variants": len(empirical),
                "mean_empirical_dosage_r2": float(empirical.mean()) if not empirical.empty else float("nan"),
                "median_empirical_dosage_r2": float(empirical.median()) if not empirical.empty else float("nan"),
                "dose0_variants": len(dose0),
                "mean_dose0": float(dose0.mean()) if not dose0.empty else float("nan"),
                "median_dose0": float(dose0.median()) if not dose0.empty else float("nan"),
                "mean_abs_dosage_bias": float(bias.mean()) if not bias.empty else float("nan"),
            }
        )
    return [{column: row.get(column, "") for column in columns} for row in rows]


def _first_existing_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    for column in candidates:
        if column in frame.columns:
            return column
    return None


def _mean_median_label(frame: pd.DataFrame, column: str | None) -> str:
    if column is None or frame.empty:
        return "not available"
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    if values.empty:
        return "not available"
    return f"{format_numeric(values.mean())} / {format_numeric(values.median())}"


def run_phasing_quality(analysis_root: Path, study: str, reference_dir: Path, force: bool) -> None:
    del reference_dir, force
    logger, _ = setup_logger(REPO_ROOT, PIPELINE_NAME, "phasing_quality")
    check_dependencies(logger, commands=[], python_packages=STAGE2_PY_PACKAGES)
    dirs = _report_dirs(analysis_root, study)
    flags = _note_flags(study, "phasing_quality")
    rows: list[dict[str, object]] = []
    assets: list[dict[str, object]] = []
    missing_chroms: list[str] = []
    proxy_chroms: list[str] = []

    for table in _stage2_phasing_metric_tables(analysis_root, study):
        frame = pd.read_csv(table, sep="\t")
        for record in frame.to_dict("records"):
            chrom = str(record.get("chrom", _chrom_from_path(table) or "NA"))
            metric_type = record.get("metric_type", "missing")
            rows.append(
                {
                    "chrom": chrom,
                    "metric_type": metric_type,
                    "variant_count": record.get("variant_count", ""),
                    "scored_variants": record.get("scored_variants", ""),
                    "mean_conf": record.get("mean_conf", ""),
                    "p5_conf": record.get("p5_conf", ""),
                    "n_low_conf_variants": record.get("n_low_conf_variants", 0),
                }
            )
            if metric_type == "ps_switch_proxy":
                proxy_chroms.append(chrom)
            elif metric_type == "missing":
                missing_chroms.append(chrom)

    if proxy_chroms:
        flags.append(
            flag_row(
                study=study,
                task="phasing_quality",
                metric="phase_set_proxy_chromosomes",
                value=",".join(sorted(set(proxy_chroms), key=_chrom_sort_key)),
                threshold="direct Eagle quality field preferred",
                flag_level="WARN",
                message="Used phase-set switch density as the phasing quality proxy for one or more chromosomes.",
            )
        )
    if missing_chroms:
        flags.append(
            flag_row(
                study=study,
                task="phasing_quality",
                metric="missing_phasing_quality_chromosomes",
                value=",".join(sorted(set(missing_chroms), key=_chrom_sort_key)),
                threshold="quality or PS field available",
                flag_level="WARN",
                message="No direct phasing quality field or PS proxy was available for one or more chromosomes.",
            )
        )

    for table in _stage2_phasing_variant_tables(analysis_root, study):
        chrom = _chrom_from_path(table) or "NA"
        assets.append(
            _asset(
                study=study,
                task="phasing_quality",
                section="Phasing",
                asset_id=f"chr{chrom}_phasing_quality_table",
                asset_type="table",
                fmt="tsv",
                path=table,
                title=f"chr{chrom} phasing quality table",
                caption="Per-variant phasing quality metric or phase-set-switch proxy.",
                sort_order=10,
            )
        )

    summary_path = dirs["tables"] / "phasing_quality_summary.tsv"
    rows = sorted(rows, key=lambda row: _chrom_sort_key(row.get("chrom", "")))
    headers = ["chrom", "metric_type", "variant_count", "scored_variants", "mean_conf", "p5_conf", "n_low_conf_variants"]
    formatted_rows = [{k: format_numeric(v) if k in ("mean_conf", "p5_conf") else v for k, v in row.items()} for row in rows]
    write_tsv(formatted_rows, summary_path, headers)
    assets.append(
        _asset(
            study=study,
            task="phasing_quality",
            section="Phasing",
            asset_id="phasing_quality_summary",
            asset_type="table",
            fmt="tsv",
            path=summary_path,
            title="Phasing quality summary",
            caption="Per-chromosome phasing quality summary derived from staged phasing metrics.",
            sort_order=11,
        )
    )
    _write_manifest(dirs, "phasing_quality", assets)
    _write_flags(dirs, "phasing_quality", flags)


def run_phasing_plots(analysis_root: Path, study: str, reference_dir: Path, force: bool) -> None:
    del reference_dir, force
    logger, _ = setup_logger(REPO_ROOT, PIPELINE_NAME, "phasing_plots")
    check_dependencies(logger, commands=[], python_packages=STAGE2_PY_PACKAGES)
    dirs = _report_dirs(analysis_root, study)
    flags = _note_flags(study, "phasing_plots")
    assets: list[dict[str, object]] = []
    frames = []

    for table in _stage2_phasing_variant_tables(analysis_root, study):
        if table.stat().st_size:
            frame = pd.read_csv(table, sep="\t")
            if not frame.empty:
                frames.append(frame)

    if not frames:
        flags.append(
            flag_row(
                study=study,
                task="phasing_plots",
                metric="phase_quality_rows",
                value=0,
                threshold=">0",
                flag_level="WARN",
                message="No per-variant phasing quality rows were available for plotting.",
            )
        )
        _write_manifest(dirs, "phasing_plots", assets)
        _write_flags(dirs, "phasing_plots", flags)
        return

    combined = pd.concat(frames, ignore_index=True)
    combined["metric_value"] = pd.to_numeric(combined.get("metric_value"), errors="coerce")
    finite = combined.dropna(subset=["metric_value"]).copy()
    if finite.empty:
        flags.append(
            flag_row(
                study=study,
                task="phasing_plots",
                metric="metric_value",
                value="all_missing",
                threshold="at least one finite value",
                flag_level="WARN",
                message="Skipped phasing plots because no finite phasing quality values were available.",
            )
        )
        _write_manifest(dirs, "phasing_plots", assets)
        _write_flags(dirs, "phasing_plots", flags)
        return

    finite["chrom_numeric"] = finite["chrom"].astype(str).str.replace("chr", "", regex=False).replace({"X": "23"}).astype(float)
    finite["global_pos"] = finite["chrom_numeric"] * 1_000_000_000 + finite["pos"].astype(float)
    fig_dir = dirs["figures"] / "phasing"

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(finite["global_pos"], finite["metric_value"], s=4, alpha=0.5)
    ax.set_title("Phasing confidence across the genome")
    ax.set_xlabel("Genome position")
    ax.set_ylabel("Confidence / proxy")
    save_figure(fig, fig_dir / "phasing_confidence_manhattan")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(finite["metric_value"], bins=50, color="steelblue", alpha=0.8)
    ax.set_title("Phasing confidence distribution")
    ax.set_xlabel("Confidence / proxy")
    ax.set_ylabel("Variants")
    save_figure(fig, fig_dir / "phasing_confidence_histogram")
    plt.close(fig)

    assets.extend(
        [
            _asset(
                study=study,
                task="phasing_plots",
                section="Phasing",
                asset_id="phasing_confidence_manhattan",
                asset_type="figure",
                fmt="png",
                path=fig_dir / "phasing_confidence_manhattan.png",
                title="Phasing confidence Manhattan plot",
                caption="Genome-wide phasing confidence or phase-set-switch proxy.",
                sort_order=20,
            ),
            _asset(
                study=study,
                task="phasing_plots",
                section="Phasing",
                asset_id="phasing_confidence_histogram",
                asset_type="figure",
                fmt="png",
                path=fig_dir / "phasing_confidence_histogram.png",
                title="Phasing confidence histogram",
                caption="Distribution of per-variant phasing confidence or proxy scores.",
                sort_order=21,
            ),
        ]
    )
    _write_manifest(dirs, "phasing_plots", assets)
    _write_flags(dirs, "phasing_plots", flags)


def run_variant_summary(analysis_root: Path, study: str, reference_dir: Path, force: bool) -> None:
    del reference_dir, force
    logger, _ = setup_logger(REPO_ROOT, PIPELINE_NAME, "variant_summary")
    check_dependencies(logger, commands=[], python_packages=STAGE2_PY_PACKAGES)
    dirs = _report_dirs(analysis_root, study)
    frame, table_path = _load_variant_frame(analysis_root, study)
    flags = _note_flags(study, "variant_summary")
    assets: list[dict[str, object]] = []
    if table_path is None or frame.empty:
        flags.append(flag_row(study=study, task="variant_summary", metric="r2_table", value="missing", threshold="present", flag_level="ERROR", message="No staged imputation R2 table was found."))
        _write_manifest(dirs, "variant_summary", assets)
        _write_flags(dirs, "variant_summary", flags)
        return

    scored = frame.dropna(subset=["r2"]).copy()
    summary_path = dirs["tables"] / "r2_summary.tsv"
    variant_summary_path = dirs["tables"] / "variant_summary.tsv"
    genotyped_maf_path = dirs["tables"] / "genotyped_maf_summary.tsv"
    imputed_maf_path = dirs["tables"] / "imputed_maf_summary.tsv"
    empirical_summary_path = dirs["tables"] / "empirical_validation_summary.tsv"
    counts, _ = np.histogram(scored["r2"].clip(0, 1), bins=R2_BINS)
    row = {
        "study": study,
        "n_total_variants": len(scored),
        "n_r2_ge_0.3": int((scored["r2"] >= 0.3).sum()),
        "n_r2_ge_0.8": int((scored["r2"] >= 0.8).sum()),
        "mean_r2": float(scored["r2"].mean()) if not scored.empty else float("nan"),
        "median_r2": float(scored["r2"].median()) if not scored.empty else float("nan"),
        "p5_r2": float(scored["r2"].quantile(0.05)) if not scored.empty else float("nan"),
    }
    row.update({label: int(value) for label, value in zip(R2_BIN_LABELS, counts, strict=True)})
    formatted_row = {k: format_numeric(v) if k in ("mean_r2", "median_r2", "p5_r2") else v for k, v in row.items()}
    write_tsv([formatted_row], summary_path)

    basis_counts = _phasing_basis_counts(analysis_root, study)
    variant_rows: list[dict[str, object]] = []
    for chrom, chrom_frame in sorted(scored.groupby("chrom_clean"), key=lambda item: _chrom_sort_key(item[0])):
        final_count = len(chrom_frame)
        classified_basis = int(chrom_frame["is_genotyped"].sum()) if "is_genotyped" in chrom_frame else 0
        basis = classified_basis or basis_counts.get(chrom, 0)
        target = max(final_count - basis, 0)
        variant_rows.append(
            {
                "chromosome": chrom,
                "imputation_basis": basis,
                "imputation_target": target,
                "total": basis + target,
            }
        )
    if variant_rows:
        variant_rows.append(
            {
                "chromosome": "Total",
                "imputation_basis": sum(int(row["imputation_basis"]) for row in variant_rows),
                "imputation_target": sum(int(row["imputation_target"]) for row in variant_rows),
                "total": sum(int(row["total"]) for row in variant_rows),
            }
        )
    write_tsv(variant_rows, variant_summary_path, ["chromosome", "imputation_basis", "imputation_target", "total"])

    genotyped = frame.loc[frame["is_genotyped"]].dropna(subset=["maf"]).copy()
    maf_rows = []
    empirical_column = _first_existing_column(genotyped, ["empirical_dosage_r2", "empirical_r2", "loo_r2", "dosage_r2"])
    dose0_column = _first_existing_column(genotyped, ["dose0", "dose_0", "mean_dose0"])
    for label, mask in [
        ("MAF < 0.05", genotyped["maf"] < 0.05 if not genotyped.empty else pd.Series(dtype=bool)),
        ("MAF >= 0.05", genotyped["maf"] >= 0.05 if not genotyped.empty else pd.Series(dtype=bool)),
    ]:
        subset = genotyped.loc[mask].copy() if not genotyped.empty else genotyped.copy()
        # Restrict count to overlapping variants with metrics
        valid_subset = subset.dropna(subset=[empirical_column, dose0_column], how="all") if empirical_column or dose0_column else subset
        maf_rows.append(
            {
                "maf_bin": label,
                "n_variants": len(valid_subset),
                "empirical_dosage_r2_mean_median": _mean_median_label(subset, empirical_column),
                "dose0_mean_median": _mean_median_label(subset, dose0_column),
            }
        )
    write_tsv(maf_rows, genotyped_maf_path, ["maf_bin", "n_variants", "empirical_dosage_r2_mean_median", "dose0_mean_median"])

    imputed = scored.loc[scored["is_imputed"]].dropna(subset=["maf"]).copy()
    imputed_maf_rows = []
    for label, mask in [
        ("MAF < 0.05", imputed["maf"] < 0.05 if not imputed.empty else pd.Series(dtype=bool)),
        ("MAF >= 0.05", imputed["maf"] >= 0.05 if not imputed.empty else pd.Series(dtype=bool)),
    ]:
        subset = imputed.loc[mask] if not imputed.empty else imputed
        imputed_maf_rows.append(
            {
                "maf_bin": label,
                "n_variants": len(subset),
                "imputation_r2_mean_median": _mean_median_label(subset, "r2"),
            }
        )
    write_tsv(imputed_maf_rows, imputed_maf_path, ["maf_bin", "n_variants", "imputation_r2_mean_median"])

    validation_mask = pd.Series(False, index=frame.index)
    for column in ("empirical_dosage_r2", "dose0", "n_samples"):
        if column in frame:
            validation_mask = validation_mask | frame[column].notna()
    empirical_frame = frame.loc[frame["is_genotyped"] & validation_mask].copy()
    empirical_rows = _empirical_validation_summary_rows(empirical_frame)
    if empirical_rows:
        formatted_rows = []
        for emp_row in empirical_rows:
            formatted_row = {
                k: format_numeric(v) if k not in ("chromosome", "empirical_r2_variants", "dose0_variants") else v
                for k, v in emp_row.items()
            }
            formatted_rows.append(formatted_row)
        write_tsv(
            formatted_rows,
            empirical_summary_path,
            [
                "chromosome",
                "empirical_r2_variants",
                "mean_empirical_dosage_r2",
                "median_empirical_dosage_r2",
                "dose0_variants",
                "mean_dose0",
                "median_dose0",
                "mean_abs_dosage_bias",
            ],
        )
    elif empirical_summary_path.exists():
        empirical_summary_path.unlink()

    pct_ge_03 = 100.0 * row["n_r2_ge_0.3"] / max(row["n_total_variants"], 1)
    if row["mean_r2"] < 0.7:
        flags.append(flag_row(study=study, task="variant_summary", metric="mean_r2", value=f"{row['mean_r2']:.4f}", threshold=">=0.7", flag_level="WARN", message="Mean imputation R2 is below 0.7."))
    if pct_ge_03 < 80.0:
        flags.append(flag_row(study=study, task="variant_summary", metric="pct_r2_ge_0.3", value=f"{pct_ge_03:.2f}", threshold=">=80%", flag_level="WARN", message="Less than 80% of variants exceed the primary R2 threshold."))
    if genotyped.empty:
        flags.append(flag_row(study=study, task="variant_summary", metric="genotyped_variant_match", value=0, threshold=">0", flag_level="WARN", message="No genotyped variants could be matched back to the Stage1 BIM; genotyped-only summaries may be unavailable."))

    assets.extend(
        [
            _asset(study=study, task="variant_summary", section="Imputation", asset_id="variant_summary", asset_type="table", fmt="tsv", path=variant_summary_path, title="Variant summary", caption="Per-chromosome counts for genotyped SNPs used as the imputation basis, SNPs generated as imputation targets, and their total.", sort_order=30),
            _asset(study=study, task="variant_summary", section="Imputation", asset_id="genotyped_maf_summary", asset_type="table", fmt="tsv", path=genotyped_maf_path, title="Genotyped variant MAF summary", caption="Genotyped validation variants grouped by MAF threshold, with empirical dosage R2 and Dose0 fields shown for overlapping sites.", sort_order=31),
            _asset(study=study, task="variant_summary", section="Imputation", asset_id="imputed_maf_summary", asset_type="table", fmt="tsv", path=imputed_maf_path, title="Imputed variant MAF summary", caption="Imputed-only variants grouped by MAF threshold, with imputation R2 shown as mean / median.", sort_order=32),
            _asset(study=study, task="variant_summary", section="Imputation", asset_id="r2_by_variant", asset_type="table", fmt="tsv", path=table_path, title="R2 by variant", caption="Per-variant imputation R2 and MAF from staged metrics across all chromosomes.", sort_order=34),
            _asset(study=study, task="variant_summary", section="Imputation", asset_id="r2_summary", asset_type="table", fmt="tsv", path=summary_path, title="R2 summary", caption="Summary metrics for imputation quality aggregated across all imputed variants.", sort_order=35),
        ]
    )
    if empirical_rows:
        assets.append(
            _asset(
                study=study,
                task="variant_summary",
                section="Imputation",
                asset_id="empirical_validation_summary",
                asset_type="table",
                fmt="tsv",
                path=empirical_summary_path,
                title="Empirical validation summary",
                caption="Genotyped validation variants compared against final imputed dosages. Empirical R2 measures dosage correlation; Dose0 is the mean imputed dosage among samples with genotype 0.",
                sort_order=33,
            )
        )
    _write_manifest(dirs, "variant_summary", assets)
    _write_flags(dirs, "variant_summary", flags)


def run_r2_by_maf(analysis_root: Path, study: str, reference_dir: Path, force: bool) -> None:
    del reference_dir, force
    logger, _ = setup_logger(REPO_ROOT, PIPELINE_NAME, "r2_by_maf")
    check_dependencies(logger, commands=[], python_packages=STAGE2_PY_PACKAGES)
    dirs = _report_dirs(analysis_root, study)
    frame, table_path = _load_variant_frame(analysis_root, study)
    flags = _note_flags(study, "r2_by_maf")
    assets: list[dict[str, object]] = []
    if table_path is None or frame.empty:
        _write_manifest(dirs, "r2_by_maf", assets)
        _write_flags(dirs, "r2_by_maf", flags)
        return

    imputed = frame.loc[frame["is_imputed"]].dropna(subset=["r2", "maf", "pos_numeric"]).copy()
    if imputed.empty:
        flags.append(flag_row(study=study, task="r2_by_maf", metric="imputed_variant_rows", value=0, threshold=">0", flag_level="WARN", message="No imputed-only variant rows were available for R2-by-MAF plots."))
        _write_manifest(dirs, "r2_by_maf", assets)
        _write_flags(dirs, "r2_by_maf", flags)
        return

    bin_summary = _maf_bin_summary(imputed, "r2")
    bin_path = dirs["tables"] / "r2_by_maf_bins.tsv"
    write_tsv(bin_summary[["maf_bin", "maf_midpoint", "n_variants", "mean_r2"]].to_dict("records"), bin_path, ["maf_bin", "maf_midpoint", "n_variants", "mean_r2"])

    fig_dir = dirs["figures"] / "imputation"
    genotyped = frame.loc[frame["is_genotyped"]].dropna(subset=["r2", "maf"]).copy()
    genotyped_summary = _maf_bin_summary(genotyped, "r2") if not genotyped.empty else pd.DataFrame()
    if genotyped.empty:
        flags.append(flag_row(study=study, task="r2_by_maf", metric="genotyped_variant_rows", value=0, threshold=">0", flag_level="WARN", message="No genotyped-only variant rows were available for the lower R2-by-MAF panel."))

    validation_mask = pd.Series(False, index=frame.index)
    for column in ("empirical_dosage_r2", "dose0", "n_samples"):
        if column in frame:
            validation_mask = validation_mask | frame[column].notna()
    validation = frame.loc[frame["is_genotyped"] & validation_mask].dropna(subset=["maf"]).copy()
    empirical_bin_path = dirs["tables"] / "empirical_validation_by_maf.tsv"

    def plot_maf_panel(ax, summary: pd.DataFrame, column: str, title: str, ylabel: str, *, r2_axis: bool = False) -> None:
        if summary.empty or column not in summary:
            ax.text(0.5, 0.5, "No rows available", ha="center", va="center", transform=ax.transAxes)
        else:
            ax.scatter(summary["maf_midpoint"], summary[column], color="black", s=20)
        ax.set_xlim(0, 0.5)
        if r2_axis:
            ax.set_ylim(0, 1.0)
        ax.set_title(title)
        ax.set_xlabel("MAF")
        ax.set_ylabel(ylabel)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    plot_maf_panel(axes[0, 0], bin_summary, "mean_r2", "A. Imputed SNP mean R2 by MAF", "Mean R2", r2_axis=True)
    plot_maf_panel(axes[0, 1], bin_summary, "n_variants", "B. Imputed SNP count by MAF", "Variant count")
    plot_maf_panel(axes[1, 0], genotyped_summary, "mean_r2", "C. Genotyped SNP mean R2 by MAF", "Mean R2", r2_axis=True)
    plot_maf_panel(axes[1, 1], genotyped_summary, "n_variants", "D. Genotyped SNP count by MAF", "Variant count")
    save_figure(fig, fig_dir / "r2_by_maf")
    plt.close(fig)

    chrom_order = sorted(imputed["chrom_clean"].dropna().unique(), key=_chrom_sort_key)
    violin_frames = []
    for chrom in chrom_order:
        chrom_frame = imputed.loc[imputed["chrom_clean"] == chrom, ["chrom_clean", "r2"]].dropna()
        if len(chrom_frame) > 50_000:
            chrom_frame = chrom_frame.sample(n=50_000, random_state=1)
        violin_frames.append(chrom_frame)
    violin_data = pd.concat(violin_frames, ignore_index=True) if violin_frames else pd.DataFrame()
    palette = {chrom: ("lightgrey" if index % 2 == 0 else "darkgrey") for index, chrom in enumerate(chrom_order)}
    fig, ax = plt.subplots(figsize=(10.5, 4.6), constrained_layout=True)
    if violin_data.empty:
        ax.text(0.5, 0.5, "No rows available", ha="center", va="center", transform=ax.transAxes)
    else:
        sns.violinplot(
            data=violin_data,
            x="chrom_clean",
            y="r2",
            hue="chrom_clean",
            order=chrom_order,
            hue_order=chrom_order,
            palette=palette,
            dodge=False,
            legend=False,
            cut=0,
            inner="quartile",
            linewidth=0.7,
            ax=ax,
        )
    ax.set_ylim(0, 1.0)
    ax.set_title("Imputed SNP R2 by chromosome")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("R2")
    ax.tick_params(axis="x", labelsize=8)
    save_figure(fig, fig_dir / "r2_by_chromosome_violin")
    plt.close(fig)

    empirical_bin_summary = _empirical_maf_bin_summary(validation)
    if empirical_bin_summary.empty:
        if empirical_bin_path.exists():
            empirical_bin_path.unlink()
        empirical_figure = fig_dir / "empirical_validation_by_maf.png"
        if empirical_figure.exists():
            empirical_figure.unlink()
    else:
        formatted_records = []
        for record in empirical_bin_summary.to_dict("records"):
            formatted_record = {
                k: format_numeric(v) if k not in ("maf_bin", "maf_midpoint", "n_variants", "n_empirical_r2", "n_dose0") else v
                for k, v in record.items()
            }
            formatted_records.append(formatted_record)

        write_tsv(
            formatted_records,
            empirical_bin_path,
            [
                "maf_bin",
                "maf_midpoint",
                "n_variants",
                "n_empirical_r2",
                "mean_empirical_dosage_r2",
                "median_empirical_dosage_r2",
                "n_dose0",
                "mean_dose0",
                "median_dose0",
            ],
        )

        def plot_empirical_panel(ax, column: str, title: str, ylabel: str, *, r2_axis: bool = False) -> None:
            data = empirical_bin_summary.dropna(subset=[column]).copy()
            if data.empty:
                ax.text(0.5, 0.5, "No rows available", ha="center", va="center", transform=ax.transAxes)
            else:
                ax.scatter(data["maf_midpoint"], data[column], color="black", s=20)
            ax.set_xlim(0, 0.5)
            if r2_axis:
                ax.set_ylim(0, 1.0)
            else:
                ax.set_ylim(bottom=0)
            ax.set_title(title)
            ax.set_xlabel("MAF")
            ax.set_ylabel(ylabel)

        fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.3), constrained_layout=True)
        plot_empirical_panel(axes[0], "mean_empirical_dosage_r2", "A. Empirical R2 by MAF", "Mean empirical R2", r2_axis=True)
        plot_empirical_panel(axes[1], "mean_dose0", "B. Dose0 by MAF", "Mean Dose0")
        plot_empirical_panel(axes[2], "n_variants", "C. Validation count by MAF", "Variant count")
        save_figure(fig, fig_dir / "empirical_validation_by_maf")
        plt.close(fig)

    assets.extend(
        [
            _asset(study=study, task="r2_by_maf", section="Imputation", asset_id="r2_by_maf_bins", asset_type="table", fmt="tsv", path=bin_path, title="Mean R2 by MAF bin", caption="Mean imputation R2 per 0.01 MAF bin for imputed-only SNPs.", sort_order=40),
            _asset(study=study, task="r2_by_maf", section="Imputation", asset_id="r2_by_maf", asset_type="figure", fmt="png", path=fig_dir / "r2_by_maf.png", title="SNP R2 by MAF", caption="Imputed SNPs (top) and genotyped SNPs (bottom) binned by 0.01 MAF intervals. Left: mean R2; Right: variant counts.", sort_order=41),
            _asset(study=study, task="r2_by_maf", section="Imputation", asset_id="r2_by_chromosome_violin", asset_type="figure", fmt="png", path=fig_dir / "r2_by_chromosome_violin.png", title="Imputed SNP R2 by chromosome", caption="Violin plots of R2 distributions by chromosome for imputed SNPs. Colors alternate between light and dark grey.", sort_order=42),
        ]
    )
    if not empirical_bin_summary.empty:
        assets.extend(
            [
                _asset(study=study, task="r2_by_maf", section="Imputation", asset_id="empirical_validation_by_maf_bins", asset_type="table", fmt="tsv", path=empirical_bin_path, title="Empirical validation by MAF bin", caption="Empirical dosage R2 and Dose0 summaries for genotyped validation variants in 0.01 MAF bins.", sort_order=43),
                _asset(study=study, task="r2_by_maf", section="Imputation", asset_id="empirical_validation_by_maf", asset_type="figure", fmt="png", path=fig_dir / "empirical_validation_by_maf.png", title="Empirical validation by MAF", caption="Panel A: empirical dosage R2 by MAF; Panel B: Dose0 by MAF; Panel C: genotyped validation variant counts.", sort_order=44),
            ]
        )
    _write_manifest(dirs, "r2_by_maf", assets)
    _write_flags(dirs, "r2_by_maf", flags)


def _reference_af_cache(analysis_root: Path, reference_dir: Path, chroms: list[str], logger) -> dict[str, Path]:
    cache_dir = ensure_dir(analysis_root / "cohort" / STAGE_NAME / "reference" / "af_cache")
    legacy_cache_dir = analysis_root / "cohort" / STAGE_NAME / "report" / "tables" / "reference_cache"
    caches: dict[str, Path] = {}
    for chrom in chroms:
        cache_path = cache_dir / f"1kg_chr{chrom}_af.tsv.gz"
        legacy_cache_path = legacy_cache_dir / cache_path.name
        if not cache_path.exists():
            if legacy_cache_path.exists():
                shutil.copy2(legacy_cache_path, cache_path)
                caches[chrom] = cache_path
                continue
            if shutil.which("bcftools") is None:
                logger.warning("bcftools is not available and no cached 1KG AF table exists for chr%s", chrom)
                continue
            vcf_path = find_reference_vcf(reference_dir, chrom)
            if vcf_path is None:
                logger.warning("Reference VCF for chr%s not found under %s", chrom, reference_dir)
                continue
            result = run_command(
                ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n", str(vcf_path)],
                logger,
                log_output=False,
            )
            with gzip.open(cache_path, "wt") as handle:
                handle.write("chrom\tpos\tref\talt\taf\n")
                handle.write(result.stdout)
        caches[chrom] = cache_path
    return caches


def run_af_concordance(analysis_root: Path, study: str, reference_dir: Path, force: bool) -> None:
    del force
    logger, _ = setup_logger(REPO_ROOT, PIPELINE_NAME, "af_concordance")
    check_dependencies(logger, commands=[], python_packages=STAGE2_PY_PACKAGES)
    dirs = _report_dirs(analysis_root, study)
    table_path = _load_or_build_r2_table(analysis_root, study)
    flags = _note_flags(study, "af_concordance")
    assets: list[dict[str, object]] = []
    if table_path is None:
        _write_manifest(dirs, "af_concordance", assets)
        _write_flags(dirs, "af_concordance", flags)
        return

    source = pd.read_csv(table_path, sep="\t").dropna(subset=["maf"])
    rows = []
    chroms = set()
    for chrom, pos, variant_id, maf in source[["chrom", "pos", "variant_id", "maf"]].itertuples(index=False, name=None):
        parts = str(variant_id).split(":")
        if len(parts) >= 4:
            chrom_key = parts[0].replace("chr", "")
            allele_key = ":".join(sorted([parts[-2].upper(), parts[-1].upper()]))
            pair_key = f"{chrom_key}:{parts[1]}:{allele_key}"
        else:
            chrom_key = str(chrom).replace("chr", "")
            pair_key = f"{chrom_key}:{pos}"
        chroms.add(chrom_key)
        rows.append({"variant_id": variant_id, "pair_key": pair_key, "imputed_af": safe_float(maf)})
    frame = pd.DataFrame(rows).dropna(subset=["imputed_af"])
    cache_paths = _reference_af_cache(analysis_root, reference_dir, sorted(chroms), logger)
    reference_rows = read_reference_af_rows(cache_paths)
    ref = pd.DataFrame(reference_rows)
    if ref.empty:
        flags.append(flag_row(study=study, task="af_concordance", metric="reference_af", value="missing", threshold="present", flag_level="WARN", message="No 1KG reference AF rows were available."))
        _write_manifest(dirs, "af_concordance", assets)
        _write_flags(dirs, "af_concordance", flags)
        return

    required_ref_columns = ["chrom", "pos", "ref", "alt", "af"]
    missing_ref_columns = [column for column in required_ref_columns if column not in ref.columns]
    if missing_ref_columns:
        flags.append(flag_row(study=study, task="af_concordance", metric="reference_af_columns", value=",".join(missing_ref_columns), threshold="chrom,pos,ref,alt,af", flag_level="WARN", message="Reference AF cache was missing required columns."))
        _write_manifest(dirs, "af_concordance", assets)
        _write_flags(dirs, "af_concordance", flags)
        return
    ref = ref.dropna(subset=[column for column in required_ref_columns if column in ref.columns]).copy()
    ref["pos_numeric"] = pd.to_numeric(ref["pos"], errors="coerce")
    ref["kg_af"] = ref["af"].astype(str).str.split(",").str[0].map(safe_float)
    ref = ref.dropna(subset=["pos_numeric", "kg_af"]).copy()
    ref["chrom_clean"] = ref["chrom"].astype(str).str.replace("chr", "", regex=False)
    ref = ref.loc[ref["chrom_clean"].ne("") & ref["ref"].astype(str).ne("") & ref["alt"].astype(str).ne("")].copy()
    ref["pair_key"] = ref.apply(lambda row: f"{row['chrom_clean']}:{int(row['pos_numeric'])}:{':'.join(sorted([str(row['ref']).upper(), str(row['alt']).upper()]))}", axis=1)
    ref["kg_af"] = ref["kg_af"].astype(float).map(lambda value: min(value, 1.0 - value))
    matched = frame.merge(ref[["pair_key", "kg_af"]], on="pair_key", how="inner")
    if matched.empty:
        flags.append(flag_row(study=study, task="af_concordance", metric="matched_variants", value=0, threshold=">0", flag_level="WARN", message="No variants matched the 1KG reference AF cache."))
        _write_manifest(dirs, "af_concordance", assets)
        _write_flags(dirs, "af_concordance", flags)
        return

    matched["af_diff"] = (matched["imputed_af"] - matched["kg_af"]).abs()
    pearson_r = float(matched[["imputed_af", "kg_af"]].corr().iloc[0, 1]) if len(matched) > 1 else float("nan")
    slope = float(np.polyfit(matched["kg_af"], matched["imputed_af"], 1)[0]) if len(matched) > 1 else float("nan")
    n_outliers = int((matched["af_diff"] > 0.10).sum())

    concordance_path = dirs["tables"] / "af_concordance.tsv"
    summary_path = dirs["tables"] / "af_concordance_summary.tsv"
    summary_path = dirs["tables"] / "af_concordance_summary.tsv"
    formatted_matched = []
    for record in matched[["variant_id", "imputed_af", "kg_af", "af_diff"]].to_dict("records"):
        formatted_record = {k: format_numeric(v) if k != "variant_id" else v for k, v in record.items()}
        formatted_matched.append(formatted_record)
    write_tsv(formatted_matched, concordance_path, ["variant_id", "imputed_af", "kg_af", "af_diff"])

    summary_row = {
        "study": study,
        "n_variants": len(matched),
        "pearson_r": pearson_r,
        "slope": slope,
        "n_outlier_variants": n_outliers,
    }
    formatted_summary = {k: format_numeric(v) if k in ("pearson_r", "slope") else v for k, v in summary_row.items()}
    write_tsv([formatted_summary], summary_path)

    fig_dir = dirs["figures"] / "imputation"
    fig, ax = plt.subplots(figsize=(4.8, 4.8), constrained_layout=True)
    scatter = ax.scatter(matched["imputed_af"], matched["kg_af"], c=matched["af_diff"], cmap="coolwarm", s=5, alpha=0.45)
    sns.regplot(data=matched, x="imputed_af", y="kg_af", scatter=False, ax=ax, color="black")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.5)
    ax.set_title(f"Allele frequency concordance\nR2={pearson_r ** 2:.4f} slope={slope:.4f}")
    ax.set_xlabel("STUDY")
    ax.set_ylabel("REF")
    fig.colorbar(scatter, ax=ax, label="|AF diff|")
    save_figure(fig, fig_dir / "af_concordance")
    plt.close(fig)

    if not math.isnan(pearson_r) and pearson_r ** 2 < 0.99:
        flags.append(flag_row(study=study, task="af_concordance", metric="r_squared", value=f"{pearson_r ** 2:.4f}", threshold=">=0.99", flag_level="WARN", message="Post-imputation AF concordance with 1KG is lower than expected."))
    assets.extend(
        [
            _asset(study=study, task="af_concordance", section="Imputation", asset_id="af_concordance", asset_type="table", fmt="tsv", path=concordance_path, title="Post-imputation AF table", caption="Per-variant post-imputation allele frequency concordance against 1000G reference.", sort_order=60),
            _asset(study=study, task="af_concordance", section="Imputation", asset_id="af_concordance_summary", asset_type="table", fmt="tsv", path=summary_path, title="Post-imputation AF summary", caption="Summary of post-imputation allele frequency agreement with 1000G reference.", sort_order=61),
            _asset(study=study, task="af_concordance", section="Imputation", asset_id="af_concordance_plot", asset_type="figure", fmt="png", path=fig_dir / "af_concordance.png", title="Post-imputation AF scatter", caption="Imputed study allele frequency compared against matched 1000G reference AF values.", sort_order=62),
        ]
    )
    _write_manifest(dirs, "af_concordance", assets)
    _write_flags(dirs, "af_concordance", flags)


TASK_DISPATCH = {
    "variant_summary": run_variant_summary,
    "r2_by_maf": run_r2_by_maf,
    "af_concordance": run_af_concordance,
}


def _read_tsv_preview(path: Path, max_rows: int = 12) -> tuple[list[str], list[dict[str, str]], int]:
    if not path.exists():
        return [], [], 0
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, str]] = []
        count = 0
        for row in reader:
            count += 1
            if len(rows) < max_rows:
                rows.append(row)
        return reader.fieldnames or [], rows, count


def _html_table(headers: list[str], rows: list[dict[str, str]]) -> str:
    if not headers:
        return "<p>No rows available.</p>"
    renamed_headers = _rename_headers(headers)
    head = "".join(f"<th>{html.escape(column)}</th>" for column in renamed_headers)
    body = []
    for row in rows:
        formatted_cells = []
        for column in headers:
            val = row.get(column, "")
            # Apply thousands separators to strings that look like integers
            if val and str(val).isdigit():
                try:
                    val = f"{int(val):,}"
                except ValueError:
                    pass
            formatted_cells.append(f"<td>{html.escape(str(val))}</td>")
        body.append("<tr>" + "".join(formatted_cells) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def _table_section(title: str, caption: str, path: Path, max_rows: int = 25) -> str:
    headers, rows, count = _read_tsv_preview(path, max_rows=max_rows)
    parts = [f"<section class=\"report-block\"><h3>{html.escape(title)}</h3>", f"<p class=\"caption\">{html.escape(caption)}</p>"]
    parts.append(_html_table(headers, rows))
    if count > len(rows):
        parts.append(f"<p class=\"note\">Showing {len(rows)} of {count} rows. Full table: <code>{html.escape(_project_relative(path))}</code></p>")
    parts.append("</section>")
    return "\n".join(parts)


def _metric_cards(dirs: dict[str, Path]) -> str:
    cards: list[tuple[str, str]] = []
    r2_summary = dirs["tables"] / "r2_summary.tsv"
    if r2_summary.exists():
        _, rows, _ = _read_tsv_preview(r2_summary, max_rows=1)
        if rows:
            row = rows[0]
            cards.extend(
                [
                    ("Variants", row.get("n_total_variants", "-")),
                    ("R2 >= 0.3", row.get("n_r2_ge_0.3", "-")),
                    ("R2 >= 0.8", row.get("n_r2_ge_0.8", "-")),
                    ("Mean R2", f"{safe_float(row.get('mean_r2')):.3f}" if safe_float(row.get("mean_r2")) is not None else "-"),
                ]
            )
    af_summary = dirs["tables"] / "af_concordance_summary.tsv"
    if af_summary.exists():
        _, rows, _ = _read_tsv_preview(af_summary, max_rows=1)
        if rows:
            row = rows[0]
            pearson = safe_float(row.get("pearson_r"))
            cards.append(("AF R2", f"{pearson ** 2:.3f}" if pearson is not None else "-"))
            cards.append(("AF outliers", row.get("n_outlier_variants", "-")))
    if not cards:
        return ""
    
    metrics_path = dirs["tables"] / f"{study}_stage2_metrics.tsv"
    metric_dict = {label: str(value).replace(",", "") for label, value in cards}
    write_tsv([metric_dict], metrics_path, [c[0] for c in cards])

    html_cards = []
    for label, value in cards:
        # Add thousands separators to whole numbers in the cards
        display_value = str(value)
        if display_value.replace(",", "").isdigit():
            try:
                display_value = f"{int(display_value.replace(',', '')):,}"
            except ValueError:
                pass
        html_cards.append(f"<div class=\"metric-card\"><span>{html.escape(label)}</span><strong>{html.escape(display_value)}</strong></div>")
    
    return f"<section><h2>Key Metrics</h2><div class=\"metric-grid\">{''.join(html_cards)}</div></section>"


def _figure_html(path: Path, title: str, caption: str) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"<figure><img src=\"data:image/png;base64,{encoded}\" alt=\"{html.escape(title)}\"><figcaption><strong>{html.escape(title)}.</strong> {html.escape(caption)}</figcaption></figure>"


def write_html_report(analysis_root: Path, study: str) -> Path:
    dirs = _report_dirs(analysis_root, study)
    html_path = dirs["report"] / "report-stage2.html"
    assets_path = dirs["report"] / "stage2_report_assets.tsv"
    flags_path = dirs["report"] / "stage2_flags.tsv"
    summary_path = analysis_root / "stage2-summary.md"

    parts = [
        "<!doctype html><html><head><meta charset=\"utf-8\">",
        "<title>Stage 2 Report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:0;background:#f7f7f7;color:#222;line-height:1.45}.page{max-width:1160px;margin:0 auto;background:#fff;padding:2rem 2.4rem 3rem}h1{font-size:1.8rem;margin:.2rem 0 .4rem}h2{font-size:1.25rem;border-bottom:2px solid #222;padding-bottom:.35rem;margin-top:1.8rem}h3{font-size:1rem;margin:.8rem 0 .2rem}.lede,.caption,.note{color:#555}.caption{margin:.15rem 0 .75rem;font-size:.92rem}table{border-collapse:collapse;width:100%;margin:.65rem 0 1rem;font-size:.86rem}th,td{border:1px solid #d8d8d8;padding:.35rem .45rem;text-align:left;vertical-align:top}th{background:#efefef}.metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:.65rem;margin:.8rem 0}.metric-card{border:1px solid #d8d8d8;background:#fafafa;padding:.7rem}.metric-card span{display:block;color:#555;font-size:.78rem;text-transform:uppercase}.metric-card strong{font-size:1.25rem}.figure-grid{display:grid;grid-template-columns:1fr;gap:1.35rem;align-items:start}figure{margin:.4rem 0 1rem}img{display:block;max-width:min(100%,980px);max-height:760px;height:auto;margin:0 auto;border:1px solid #ddd}figcaption{font-size:.9rem;color:#444;margin:.4rem auto 0;max-width:980px}.report-block{margin-bottom:1.2rem}code{font-size:.86rem}.flag-warn{color:#8a5a00}.flag-error{color:#a40000}@media(max-width:760px){.page{padding:1rem}img{max-height:none}}</style>",
        "</head><body>",
        "<main class=\"page\">",
        f"<h1>Stage 2 Report: {html.escape(study)}</h1>",
        "<p class=\"lede\">This report summarizes Stage 2 phasing and imputation outputs for this study.</p>",
    ]

    if summary_path.exists():
        parts.append(f"<p class=\"note\">Authoritative cross-study summary: <code>{html.escape(_project_relative(summary_path))}</code></p>")

    metric_cards = _metric_cards(dirs)
    if metric_cards:
        parts.append(metric_cards)

    parts.append("<h2>Tables</h2>")
    variant_summary = dirs["tables"] / "variant_summary.tsv"
    if variant_summary.exists():
        parts.append(_table_section("Variant summary", "Counts of SNPs used as the imputation basis, SNPs generated as imputation targets, and their total by chromosome.", variant_summary, max_rows=30))
    genotyped_maf_summary = dirs["tables"] / "genotyped_maf_summary.tsv"
    if genotyped_maf_summary.exists():
        parts.append(_table_section("MAF of genotyped variants", "Genotyped validation variants split at MAF 0.05. Empirical dosage R2 and Dose0 are shown for overlapping sites.", genotyped_maf_summary, max_rows=5))
    imputed_maf_summary = dirs["tables"] / "imputed_maf_summary.tsv"
    if imputed_maf_summary.exists():
        parts.append(_table_section("MAF of imputed variants", "Imputed-only variants split at MAF 0.05. Imputation R2 is shown as mean / median.", imputed_maf_summary, max_rows=5))
    empirical_summary = dirs["tables"] / "empirical_validation_summary.tsv"
    if empirical_summary.exists():
        parts.append(_table_section("Empirical validation summary", "Genotyped variants compared against final imputed dosages by chromosome. Empirical R2 measures dosage correlation; Dose0 is the mean imputed dosage among samples observed with genotype 0.", empirical_summary, max_rows=30))
    if flags_path.exists():
        parts.append(_table_section("Flags", "Warnings and notes raised while preparing the Stage 2 report.", flags_path, max_rows=200))
    figure_specs = [
        (dirs["figures"] / "imputation" / "r2_by_maf.png", "SNP R2 by MAF", "Top row summarizes imputed SNPs and bottom row summarizes genotyped SNPs in 0.01 MAF bins. Left panels show mean R2 and right panels show variant counts."),
        (dirs["figures"] / "imputation" / "empirical_validation_by_maf.png", "Empirical validation by MAF", "Typed validation variants are binned by MAF. Panels show empirical dosage R2, Dose0, and variant counts."),
        (dirs["figures"] / "imputation" / "r2_by_chromosome_violin.png", "Imputed SNP R2 by chromosome", "Violin plot of R2 distributions for imputed SNPs by chromosome."),
        (dirs["figures"] / "imputation" / "af_concordance.png", "Allele frequency concordance", "Study allele frequency is on the X axis and reference allele frequency is on the Y axis."),
    ]
    figures = [(path, title, caption) for path, title, caption in figure_specs if path.exists()]
    if figures:
        parts.append("<h2>Figures</h2>")
        parts.append("<div class=\"figure-grid\">")
        for figure, title, caption in figures:
            parts.append(_figure_html(figure, title, caption))
        parts.append("</div>")

    if assets_path.exists():
        parts.append(f"<p class=\"note\">Report asset manifest: <code>{html.escape(_project_relative(assets_path))}</code></p>")

    parts.append("</main></body></html>")
    html_path.write_text("\n".join(parts))
    return html_path


def run_task(task_id: str, analysis_root: Path, studies_arg: str, reference_dir: Path, force: bool = False, chromosomes_arg: str = "all") -> None:
    global ACTIVE_CHROMS
    ACTIVE_CHROMS = _parse_chromosomes(chromosomes_arg)
    studies = _study_ids(studies_arg, analysis_root)
    for study in studies:
        _prepare_study_report_dir(analysis_root, study)
        TASK_DISPATCH[task_id](analysis_root, study, reference_dir, force)
        merge_stage_study_assets(analysis_root, STAGE_NAME, study)
        merge_stage_study_flags(analysis_root, STAGE_NAME, study)
        write_html_report(analysis_root, study)


def run_reports(analysis_root: Path, studies_arg: str, reference_dir: Path, force: bool = False, chromosomes_arg: str = "all") -> None:
    global ACTIVE_CHROMS
    ACTIVE_CHROMS = _parse_chromosomes(chromosomes_arg)
    studies = _study_ids(studies_arg, analysis_root)
    _migrate_legacy_reference_cache(analysis_root)
    cohort_report = analysis_root / "cohort" / STAGE_NAME / "report"
    if cohort_report.exists():
        shutil.rmtree(cohort_report)
    for study in studies:
        _prepare_study_report_dir(analysis_root, study)
        for task_id in TASK_SEQUENCE:
            TASK_DISPATCH[task_id](analysis_root, study, reference_dir, force)
        merge_stage_study_assets(analysis_root, STAGE_NAME, study)
        merge_stage_study_flags(analysis_root, STAGE_NAME, study)
        write_html_report(analysis_root, study)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create semantic per-study Stage 2 report artifacts.")
    parser.add_argument("task_id", nargs="?", choices=sorted(TASK_DISPATCH))
    parser.add_argument("--analysis-root", "--analysis_root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--studies", "--study", default="all")
    parser.add_argument("--chromosomes", "--chromosome", "--chrom", default="all")
    parser.add_argument("--reference-dir", "--reference_dir", default=str(REPO_ROOT / "data" / "reference" / "1000G"))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    reference_dir = Path(args.reference_dir).resolve()
    if args.task_id:
        run_task(args.task_id, analysis_root, args.studies, reference_dir, args.force, args.chromosomes)
    else:
        run_reports(analysis_root, args.studies, reference_dir, args.force, args.chromosomes)


if __name__ == "__main__":
    main()
