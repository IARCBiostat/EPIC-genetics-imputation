from __future__ import annotations

import csv
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Iterable


def q(value: str | Path) -> str:
    return shlex.quote(str(value))


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def run(cmd: str, cwd: Path) -> None:
    print(cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True, cwd=str(cwd))


def configure_plink(plink: str) -> str:
    resolved = Path(plink).expanduser() if os.path.sep in plink else None
    if resolved is not None:
        if not resolved.exists():
            raise FileNotFoundError(f"PLINK executable not found: {resolved}")
        if not os.access(resolved, os.X_OK):
            raise PermissionError(f"PLINK executable is not executable: {resolved}")
        plink_path = str(resolved.resolve())
    else:
        found = shutil.which(plink)
        if found is None:
            raise FileNotFoundError(
                f"PLINK executable not found in PATH: {plink}. "
                "Pass --plink /path/to/plink or add PLINK to PATH."
            )
        plink_path = found

    os.environ["PLINK_BIN"] = plink_path
    plink_dir = str(Path(plink_path).parent)
    path_parts = os.environ.get("PATH", "").split(os.pathsep)
    if plink_dir not in path_parts:
        os.environ["PATH"] = plink_dir + os.pathsep + os.environ.get("PATH", "")
    return plink_path


def copy_prefix(src_prefix: Path, dst_prefix: Path, cwd: Path) -> None:
    for ext in ("bim", "fam", "bed"):
        run(f"cp {q(str(src_prefix) + '.' + ext)} {q(str(dst_prefix) + '.' + ext)}", cwd)


def move_prefix(src_prefix: Path, dst_prefix: Path, cwd: Path) -> None:
    for ext in ("bim", "fam", "bed"):
        run(f"mv {q(str(src_prefix) + '.' + ext)} {q(str(dst_prefix) + '.' + ext)}", cwd)


def remove_prefix(prefix: Path, cwd: Path) -> None:
    run(f"rm -f {q(str(prefix))}.*", cwd)


def write_empty_file(path: Path) -> None:
    path.write_text("")


def write_summary(path: Path, text: str, mode: str = "a") -> None:
    with path.open(mode) as handle:
        handle.write(text)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    with path.open("w") as handle:
        handle.writelines(lines)


STUDY_METADATA_COLUMNS = {
    "Brea_01_Erneg": "Brea_01_Erneg",
    "Brea_02": "Brea_02_Onco",
    "Brea_02_Onco": "Brea_02_Onco",
    "Clrt_01": "Clrt_01_Gecco",
    "Clrt_01_Gecco": "Clrt_01_Gecco",
    "Ecvd_01": "Ecvd_01",
    "Ecvd_02": "Ecvd_02",
    "Ecvd_03": "Ecvd_03",
    "Glbd_01": "Glbd_01",
    "Inte_01": "Inte_01",
    "Inte_02": "Inte_02",
    "Inte_03": "Inte_03",
    "Kidn_01": "Kidn_01",
    "Kidn_02": "Kidn_02",
    "Lung_01": "Lung_01",
    "Lymp_01": "Lymp_01",
    "Neuro_01": "Neuro_01",
    "Ovar_01": "Ovar_01",
    "Panc_01": "Panc_01_PS1",
    "Panc_01_PS1": "Panc_01_PS1",
    "Panc_02": "Panc_02_PS3",
    "Panc_02_PS3": "Panc_02_PS3",
    "Pros_01": "Pros_01_Bpc3",
    "Pros_01_Bpc3": "Pros_01_Bpc3",

    "Pros_03": "Pros_03_Onco",
    "Pros_03_Onco": "Pros_03_Onco",
    "Pros_04": "Pros_04_P160555",
    "Pros_04_P160555": "Pros_04_P160555",
    "Stom_01": "Stom_01",
    "Uadt_01": "Uadt_01",
}


_EPIC_SAMPLE_SUFFIX_RE = re.compile(r"(_R\d+|_QC\d*)$")


def _base_epic_sample_id(sample_id: str) -> str:
    return _EPIC_SAMPLE_SUFFIX_RE.sub("", sample_id)


def resolve_epic_reference_file(data_root: Path, filename: str) -> Path:
    if filename == "EPIC_study_case_status.txt":
        override = os.environ.get("EPIC_CASE_STATUS_FILE")
        if override:
            override_path = Path(override)
            if override_path.exists():
                return override_path
            raise FileNotFoundError(f"EPIC_CASE_STATUS_FILE is set but does not exist: {override_path}")

    candidates = [
        data_root / "reference" / "Epic" / filename,
        data_root / "Reference" / "Epic" / filename,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not find EPIC reference file "
        f"{filename!r}; tried: {', '.join(str(path) for path in candidates)}"
    )


def _load_epic_sample_metadata(metadata_path: Path, study_column: str) -> dict[str, tuple[str, str]]:
    metadata: dict[str, tuple[str, str]] = {}
    with metadata_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"EPIC metadata file is empty: {metadata_path}")
        required = {"ID", "sex", study_column}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"EPIC metadata file {metadata_path} is missing required column(s): "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            sample_id = row["ID"].strip()
            if not sample_id:
                continue
            sex = row["sex"].strip()
            pheno = row[study_column].strip()
            sex = sex if sex in {"1", "2"} else "0"
            pheno = pheno if pheno in {"1", "2"} else "-9"

            if sample_id in metadata and metadata[sample_id] != (sex, pheno):
                raise ValueError(f"Conflicting EPIC metadata rows for sample ID {sample_id}")
            metadata[sample_id] = (sex, pheno)

    return metadata


def apply_epic_sample_metadata(
    source_prefix: Path,
    output_prefix: Path,
    data_root: Path,
    study_id: str,
    summary: Path,
    strict: bool = True,
) -> None:
    study_column = STUDY_METADATA_COLUMNS.get(study_id, study_id)
    metadata_path = resolve_epic_reference_file(data_root, "EPIC_study_case_status.txt")
    metadata = _load_epic_sample_metadata(metadata_path, study_column)

    source_fam = source_prefix.with_suffix(".fam")
    output_fam = output_prefix.with_suffix(".fam")
    output_fam.parent.mkdir(parents=True, exist_ok=True)

    sample_count = 0
    sex_counts = {"0": 0, "1": 0, "2": 0}
    pheno_counts = {"-9": 0, "1": 0, "2": 0}
    unmatched: list[str] = []
    unmatched_rows: list[tuple[str, str, str]] = []

    with source_fam.open() as src, output_fam.open("w") as dst:
        for line_number, line in enumerate(src, start=1):
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"Malformed FAM row in {source_fam}:{line_number}")

            lookup_id = _base_epic_sample_id(fields[1])
            if lookup_id not in metadata:
                unmatched.append(fields[1])
                unmatched_rows.append((fields[0], fields[1], lookup_id))
                sex = "0"
                pheno = "-9"
            else:
                sex, pheno = metadata[lookup_id]

            fields[4] = sex
            fields[5] = pheno
            dst.write("\t".join(fields[:6]) + "\n")

            sample_count += 1
            sex_counts[sex] = sex_counts.get(sex, 0) + 1
            pheno_counts[pheno] = pheno_counts.get(pheno, 0) + 1

    if unmatched and strict:
        preview = ", ".join(unmatched[:10])
        suffix = " ..." if len(unmatched) > 10 else ""
        output_fam.unlink(missing_ok=True)
        raise ValueError(
            f"{len(unmatched)} sample(s) in {source_fam} were not found in {metadata_path} "
            f"using FAM IID lookup for study {study_id}: {preview}{suffix}"
        )

    for ext in ("bed", "bim"):
        shutil.copyfile(source_prefix.with_suffix(f".{ext}"), output_prefix.with_suffix(f".{ext}"))

    write_summary(
        summary,
        (
            "\n*************************** EPIC Sample Metadata *******************************\n\n"
            f"Metadata file: {metadata_path}\n"
            f"Study phenotype column: {study_column}\n"
            f"Samples updated: {sample_count}\n"
            f"Sex counts: unknown={sex_counts.get('0', 0)}, male={sex_counts.get('1', 0)}, female={sex_counts.get('2', 0)}\n"
            f"Phenotype counts: missing={pheno_counts.get('-9', 0)}, control={pheno_counts.get('1', 0)}, case={pheno_counts.get('2', 0)}\n"
            f"Unmatched samples: {len(unmatched)}\n\n"
            "**********************************************************\n\n"
        ),
        "a",
    )

    if unmatched_rows:
        unmatched_path = output_prefix.with_name(f"{output_prefix.name}_unmatched_epic_metadata.tsv")
        with unmatched_path.open("w") as handle:
            handle.write("FID\tIID\tlookup_id\n")
            for fid, iid, lookup_id in unmatched_rows:
                handle.write(f"{fid}\t{iid}\t{lookup_id}\n")
        write_summary(
            summary,
            (
                f"Unmatched sample audit file: {unmatched_path}\n"
                "Unmatched samples were retained with PLINK-missing sex=0 and phenotype=-9.\n\n"
            ),
            "a",
        )


def manifest_arg(data_root: Path, rel_path: str, flags: str) -> str:
    return f"{q(data_root / rel_path)} {flags}".strip()


def python2_cmd(python2: str, script_path: Path, args: str) -> str:
    return f"{q(python2)} {q(script_path)} {args}".strip()


def plink_cmd(plink: str, args: str) -> str:
    return f"{q(plink)} {args}".strip()


def _source_build_name(build: str) -> str:
    if build == "36":
        return "hg18"
    if build == "37":
        return "hg19"
    if build == "38":
        return "hg38"
    raise ValueError(f"Unsupported genome build for liftover: {build}")


def _normalize_target_chr(value: str) -> str:
    mapping = {
        "M": "MT",
        "MT": "MT",
        "X": "X",
        "Y": "Y",
    }
    return mapping.get(value, value)


def _prepare_liftover_files(source_prefix: Path, invr_path: Path, out_prefix: Path) -> tuple[Path, Path, Path, Path, Path, int, int]:
    original_ids: list[str] = []
    bim_alleles: dict[str, tuple[str, str]] = {}
    with source_prefix.with_suffix(".bim").open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            snpid = fields[1]
            original_ids.append(snpid)
            bim_alleles[snpid] = (fields[4], fields[5])

    lifted_ids: list[str] = []
    lifted_id_set: set[str] = set()
    chr_updates: list[str] = []
    pos_updates: list[str] = []
    inverted_snvs: list[str] = []

    with invr_path.open() as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue
            snpid, chrom, target_pos, _source_pos, category = fields[:5]
            lifted_ids.append(snpid)
            lifted_id_set.add(snpid)
            chr_updates.append(f"{snpid}\t{_normalize_target_chr(chrom)}\n")
            pos_updates.append(f"{snpid}\t{target_pos}\n")
            a1, a2 = bim_alleles.get(snpid, ("", ""))
            if category == "inverted" and len(a1) == 1 and len(a2) == 1 and a1 != "0" and a2 != "0":
                inverted_snvs.append(f"{snpid}\n")

    unlifted_ids = [f"{snpid}\n" for snpid in original_ids if snpid not in lifted_id_set]

    extract_path = out_prefix.with_name(f"{out_prefix.name}_lifted_variants.txt")
    chr_path = out_prefix.with_name(f"{out_prefix.name}_lifted_chr.txt")
    pos_path = out_prefix.with_name(f"{out_prefix.name}_lifted_pos.txt")
    flip_path = out_prefix.with_name(f"{out_prefix.name}_inverted_snvs.txt")
    unlifted_path = out_prefix.with_name(f"{out_prefix.name}_unlifted.txt")

    write_lines(extract_path, [f"{snpid}\n" for snpid in lifted_ids])
    write_lines(chr_path, chr_updates)
    write_lines(pos_path, pos_updates)
    write_lines(flip_path, inverted_snvs)
    write_lines(unlifted_path, unlifted_ids)

    return (
        extract_path,
        chr_path,
        pos_path,
        flip_path,
        unlifted_path,
        len(original_ids),
        len(lifted_ids),
    )


def liftover_to_hg38(
    plink: str,
    repo_root: Path,
    source_prefix: Path,
    output_prefix: Path,
    build: str,
    lift_dir: Path,
    summary: Path,
    cwd: Path,
) -> None:
    source_build = _source_build_name(build)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    if source_build == "hg38":
        move_prefix(source_prefix, output_prefix, cwd)
        return

    triple_dir = Path(os.environ.get("STAGE1_TRIPLE_LIFTOVER_DIR", repo_root / "tools" / "triple-liftOver")).resolve()
    triple_script = triple_dir / "tripleliftover_v133.pl"
    chain_gz = triple_dir / "library" / "chainfiles" / f"{source_build}ToHg38.over.chain.gz"
    if not triple_script.exists():
        raise FileNotFoundError(
            f"Missing triple-liftOver script: {triple_script}. "
            "Set STAGE1_TRIPLE_LIFTOVER_DIR to the server triple-liftOver installation."
        )
    if not chain_gz.exists():
        raise FileNotFoundError(
            f"Missing liftover chain file: {chain_gz}. "
            "Set STAGE1_TRIPLE_LIFTOVER_DIR to the server triple-liftOver installation."
        )
    triple_prefix = lift_dir / output_prefix.name
    run(
        f"perl {q(triple_script)} --bim {q(source_prefix.with_suffix('.bim'))} --base {source_build} --target hg38 --outprefix {q(triple_prefix)}",
        cwd,
    )

    invr_path = triple_prefix.with_suffix(".invr.txt")
    (
        extract_path,
        chr_path,
        pos_path,
        flip_path,
        _unlifted_path,
        input_count,
        lifted_count,
    ) = _prepare_liftover_files(source_prefix, invr_path, triple_prefix)

    if lifted_count == 0:
        raise RuntimeError(f"Liftover produced no variants for {source_prefix}")

    lifted_prefix = lift_dir / f"{output_prefix.name}_lifted"
    chr_prefix = lift_dir / f"{output_prefix.name}_lifted_chr"
    pos_prefix = lift_dir / f"{output_prefix.name}_lifted_chr_pos"
    flip_prefix = lift_dir / f"{output_prefix.name}_lifted_chr_pos_flip"

    run(
        plink_cmd(plink, f"--bfile {q(source_prefix)} --extract {q(extract_path)} --make-bed --out {q(lifted_prefix)}"),
        cwd,
    )
    run(
        plink_cmd(plink, f"--bfile {q(lifted_prefix)} --update-chr {q(chr_path)} 2 1 --make-bed --out {q(chr_prefix)}"),
        cwd,
    )
    remove_prefix(lifted_prefix, cwd)
    run(
        plink_cmd(plink, f"--bfile {q(chr_prefix)} --update-map {q(pos_path)} 2 1 --make-bed --out {q(pos_prefix)}"),
        cwd,
    )
    remove_prefix(chr_prefix, cwd)

    split_source = pos_prefix
    if flip_path.stat().st_size > 0:
        run(
            plink_cmd(plink, f"--bfile {q(pos_prefix)} --flip {q(flip_path)} --make-bed --out {q(flip_prefix)}"),
            cwd,
        )
        remove_prefix(pos_prefix, cwd)
        split_source = flip_prefix

    tmp_output = lift_dir / output_prefix.name
    run(
        plink_cmd(plink, f"--bfile {q(split_source)} --split-x b38 no-fail --make-bed --out {q(tmp_output)}"),
        cwd,
    )
    remove_prefix(split_source, cwd)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    move_prefix(tmp_output, output_prefix, cwd)

    write_summary(
        summary,
        (
            "\n**************************** LiftOver results ******************************\n\n"
            f"Input build: {source_build}\n"
            "Target build: hg38\n"
            f"Number of SNPs before using LiftOver: {input_count} SNPs.\n"
            f"Number of SNPs after using LiftOver: {lifted_count} SNPs.\n"
            f"Number of SNPs lost: {input_count - lifted_count} SNPs.\n\n"
            "**********************************************************\n\n"
        ),
        "a",
    )
