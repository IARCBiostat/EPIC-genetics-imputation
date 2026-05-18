from __future__ import annotations

import csv
import gzip
import importlib
import logging
import math
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")


DEFAULT_STUDIES = [
    "Brea_01_Erneg",
    "Brea_02",
    "Clrt_01",
    "Ecvd_01",
    "Ecvd_02",
    "Ecvd_03",
    "Glbd_01",
    "Inte_01",
    "Inte_02",
    "Inte_03",
    "Kidn_01",
    "Kidn_02",
    "Lung_01",
    "Lymp_01",
    "Ovar_01",
    "Panc_01",
    "Panc_02",
    "Pros_01",

    "Pros_03",
    "Pros_04",
    "Stom_01",
    "Uadt_01",
]

SPECIAL_STUDY_NOTES = {
    "Clrt_01": "Missing entirely; skip gracefully in all tasks.",
    "Kidn_02": "Variant-count outlier (~4.1M variants); flag in per-variant metrics.",
    "Glbd_01": "Small study (N=119); flag in per-sample metrics.",
    "Uadt_01": "Small study (N=213); flag in per-sample metrics.",
    "Panc_02": "Small study (N=183); flag in per-sample metrics and build-35 exclusion metadata.",
    "Inte_01": "Build 36 (hg18) origin; watch for coordinate anomalies.",
    "Kidn_01": "Build 36 (hg18) origin; watch for coordinate anomalies.",
    "Lymp_01": "Build-35 exclusion applied in Stage 1; note in QC metadata.",
}

ASSET_COLUMNS = [
    "study",
    "stage",
    "task",
    "section",
    "asset_id",
    "asset_type",
    "format",
    "path",
    "title",
    "caption",
    "sort_order",
]

FLAG_COLUMNS = [
    "study",
    "task",
    "metric",
    "value",
    "threshold",
    "flag_level",
    "message",
]


def repo_root_from(current_file: str | Path) -> Path:
    return Path(current_file).resolve().parents[2]


def parse_studies(studies_arg: str | None, analysis_root: Path | None = None) -> list[str]:
    if studies_arg and studies_arg != "all":
        return sorted({item.strip() for item in studies_arg.split(",") if item.strip()})
    if analysis_root is None:
        return DEFAULT_STUDIES.copy()

    discovered: set[str] = set()
    for stage_name in ("stage1", "stage2", "stage3"):
        for path in analysis_root.glob(f"*/{stage_name}"):
            if path.is_dir() and path.parent.name != "cohort":
                discovered.add(path.parent.name)
    return sorted(discovered or set(DEFAULT_STUDIES))


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def stage_dir(analysis_root: Path, study: str, stage_name: str) -> Path:
    return analysis_root / study / stage_name


def study_qc_dirs(analysis_root: Path, study: str, stage_name: str) -> dict[str, Path]:
    root = stage_dir(analysis_root, study, stage_name)
    report = root / "report"
    tables = report / "tables"
    figures = report / "figures"
    manifests = report / "manifests"
    flags = report / "flags"
    report_inputs = analysis_root / study / "report_inputs"
    for path in (tables, figures, manifests, flags, report_inputs):
        ensure_dir(path)
    return {
        "root": root,
        "report": report,
        "qc": report,
        "tables": tables,
        "figures": figures,
        "manifests": manifests,
        "flags": flags,
        "report_inputs": report_inputs,
    }


def cohort_qc_dirs(analysis_root: Path, stage_name: str) -> dict[str, Path]:
    root = analysis_root / "cohort" / stage_name / "report"
    tables = root / "tables"
    figures = root / "figures"
    manifests = root / "manifests"
    flags = root / "flags"
    for path in (tables, figures, manifests, flags):
        ensure_dir(path)
    return {
        "root": root,
        "tables": tables,
        "figures": figures,
        "manifests": manifests,
        "flags": flags,
    }


def task_done_path(analysis_root: Path, stage_name: str, task_id: str) -> Path:
    return ensure_dir(analysis_root / "cohort" / stage_name / "report") / f".{task_id}.done"


def write_done(path: Path, message: str = "") -> None:
    path.write_text(f"{utc_timestamp()}\n{message}".rstrip() + "\n")


def setup_logger(repo_root: Path, pipeline_name: str, task_id: str) -> tuple[logging.Logger, Path]:
    log_dir = ensure_dir(repo_root / pipeline_name / "logs")
    log_path = log_dir / f"{task_id}_{utc_timestamp()}.log"
    logger_name = f"{pipeline_name}.{task_id}.{utc_timestamp()}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger, log_path


def check_dependencies(
    logger: logging.Logger,
    commands: Sequence[str] | None = None,
    python_packages: Sequence[str] | None = None,
) -> None:
    missing_commands = [cmd for cmd in (commands or []) if shutil.which(cmd) is None]
    missing_packages = []
    for package in python_packages or []:
        try:
            importlib.import_module(package)
        except ImportError:
            missing_packages.append(package)

    if missing_commands or missing_packages:
        if missing_commands:
            logger.error("Missing required commands: %s", ", ".join(missing_commands))
        if missing_packages:
            logger.error("Missing required Python packages: %s", ", ".join(missing_packages))
        raise RuntimeError("Dependency check failed")


def run_command(
    command: Sequence[str],
    logger: logging.Logger,
    cwd: Path | None = None,
    allow_failure: bool = False,
    log_output: bool = True,
) -> subprocess.CompletedProcess[str]:
    logger.info("Running command: %s", " ".join(command))
    process = subprocess.Popen(
        list(command),
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert process.stdout is not None
    captured: list[str] = []
    for line in process.stdout:
        captured.append(line)
        if log_output:
            logger.info(line.rstrip("\n"))
    return_code = process.wait()
    output = "".join(captured)
    if not log_output:
        logger.info("Command output logging suppressed; captured %s lines.", len(captured))
    result = subprocess.CompletedProcess(command, return_code, output, "")
    if return_code != 0 and not allow_failure:
        raise subprocess.CalledProcessError(return_code, command, output=output)
    return result


def count_lines(path: Path, skip_blank: bool = True, skip_comments: bool = False) -> int:
    if not path.exists():
        return 0
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if skip_blank and not line.strip():
                continue
            if skip_comments and line.startswith("#"):
                continue
            count += 1
    return count


def write_tsv(rows: Sequence[dict[str, object]], path: Path, fieldnames: Sequence[str] | None = None) -> None:
    ensure_dir(path.parent)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge_tsv_files(input_paths: Iterable[Path], output_path: Path, columns: Sequence[str]) -> None:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    for path in input_paths:
        for row in read_tsv(path):
            key = tuple(row.get(column, "") for column in columns)
            if key in seen:
                continue
            seen.add(key)
            rows.append({column: row.get(column, "") for column in columns})
    write_tsv(rows, output_path, fieldnames=columns)


def save_figure(fig, base_path: Path, dpi: int = 300) -> Path:
    ensure_dir(base_path.parent)
    png_path = base_path.with_suffix(".png")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    return png_path


def asset_row(
    *,
    study: str,
    stage: str,
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
    return {
        "study": study,
        "stage": stage,
        "task": task,
        "section": section,
        "asset_id": asset_id,
        "asset_type": asset_type,
        "format": fmt,
        "path": str(path),
        "title": title,
        "caption": caption,
        "sort_order": sort_order,
    }


def flag_row(
    *,
    study: str,
    task: str,
    metric: str,
    value: object,
    threshold: str,
    flag_level: str,
    message: str,
) -> dict[str, object]:
    return {
        "study": study,
        "task": task,
        "metric": metric,
        "value": value,
        "threshold": threshold,
        "flag_level": flag_level,
        "message": message,
    }


def manifest_output_path(dir_map: dict[str, Path], task_id: str) -> Path:
    return dir_map["manifests"] / f"{task_id}.assets.tsv"


def flags_output_path(dir_map: dict[str, Path], task_id: str) -> Path:
    return dir_map["flags"] / f"{task_id}.flags.tsv"


def merge_stage_study_assets(analysis_root: Path, stage_name: str, study: str) -> Path:
    dir_map = study_qc_dirs(analysis_root, study, stage_name)
    output = dir_map["report"] / f"{stage_name}_report_assets.tsv"
    merge_tsv_files(sorted(dir_map["manifests"].glob("*.assets.tsv")), output, ASSET_COLUMNS)
    return output


def merge_stage_study_flags(analysis_root: Path, stage_name: str, study: str) -> Path:
    dir_map = study_qc_dirs(analysis_root, study, stage_name)
    output = dir_map["report"] / f"{stage_name}_flags.tsv"
    merge_tsv_files(sorted(dir_map["flags"].glob("*.flags.tsv")), output, FLAG_COLUMNS)
    return output


def merge_stage_cohort_assets(analysis_root: Path, stage_name: str) -> Path:
    dir_map = cohort_qc_dirs(analysis_root, stage_name)
    output = dir_map["root"] / f"{stage_name}_cohort_report_assets.tsv"
    merge_tsv_files(sorted(dir_map["manifests"].glob("*.assets.tsv")), output, ASSET_COLUMNS)
    return output


def merge_stage_cohort_flags(analysis_root: Path, stage_name: str) -> Path:
    dir_map = cohort_qc_dirs(analysis_root, stage_name)
    output = dir_map["root"] / f"{stage_name}_cohort_flags.tsv"
    merge_tsv_files(sorted(dir_map["flags"].glob("*.flags.tsv")), output, FLAG_COLUMNS)
    return output


def merge_all_stage_assets(analysis_root: Path, study: str) -> Path:
    report_dir = ensure_dir(analysis_root / study / "report_inputs")
    output = report_dir / f"{study}_all_stages_report_assets.tsv"
    input_paths = [
        analysis_root / study / "stage1" / "report" / "stage1_report_assets.tsv",
        analysis_root / study / "stage2" / "report" / "stage2_report_assets.tsv",
        analysis_root / study / "stage3" / "report" / "stage3_report_assets.tsv",
    ]
    merge_tsv_files([path for path in input_paths if path.exists()], output, ASSET_COLUMNS)
    return output


def merge_all_stage_flags(analysis_root: Path, study: str) -> Path:
    report_dir = ensure_dir(analysis_root / study / "report_inputs")
    output = report_dir / f"{study}_warnings.tsv"
    input_paths = [
        analysis_root / study / "stage1" / "report" / "stage1_flags.tsv",
        analysis_root / study / "stage2" / "report" / "stage2_flags.tsv",
        analysis_root / study / "stage3" / "report" / "stage3_flags.tsv",
    ]
    merge_tsv_files([path for path in input_paths if path.exists()], output, FLAG_COLUMNS)
    return output


def find_reference_vcf(reference_dir: Path, chrom: str) -> Path | None:
    candidates = []
    patterns = [
        f"*chr{chrom}*.vcf.gz",
        f"*CHR{chrom}*.vcf.gz",
        f"*{chrom}*.vcf.gz",
    ]
    for pattern in patterns:
        candidates.extend(sorted(reference_dir.glob(pattern)))
    for path in candidates:
        if path.name.endswith(".tbi"):
            continue
        return path
    return None


def reference_af_cache(
    analysis_root: Path,
    stage_name: str,
    reference_dir: Path,
    chroms: Sequence[str],
    logger: logging.Logger,
) -> dict[str, Path]:
    cache_dir = ensure_dir(cohort_qc_dirs(analysis_root, stage_name)["tables"] / "reference_cache")
    caches: dict[str, Path] = {}
    for chrom in chroms:
        cache_path = cache_dir / f"1kg_chr{chrom}_af.tsv.gz"
        if not cache_path.exists():
            vcf_path = find_reference_vcf(reference_dir, chrom)
            if vcf_path is None:
                logger.warning("Reference VCF for chr%s not found under %s", chrom, reference_dir)
                continue
            command = [
                shutil.which("bcftools") or "bcftools",
                "query",
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
                str(vcf_path),
            ]
            logger.info("Caching 1KG AF for chr%s from %s", chrom, vcf_path)
            result = run_command(command, logger, log_output=False)
            with gzip.open(cache_path, "wt") as handle:
                handle.write("chrom\tpos\tref\talt\taf\n")
                handle.write(result.stdout)
        caches[chrom] = cache_path
    return caches


def read_reference_af_rows(cache_paths: dict[str, Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in cache_paths.values():
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows.extend(reader)
    return rows


def safe_float(value: str | None) -> float | None:
    if value in {None, "", ".", "NA", "nan"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def is_strand_ambiguous(allele1: str, allele2: str) -> bool:
    pair = {allele1.upper(), allele2.upper()}
    return pair in ({"A", "T"}, {"C", "G"})


def cohort_median(values: Sequence[float]) -> float | None:
    filtered = sorted(value for value in values if not math.isnan(value))
    if not filtered:
        return None
    mid = len(filtered) // 2
    if len(filtered) % 2:
        return filtered[mid]
    return (filtered[mid - 1] + filtered[mid]) / 2.0
