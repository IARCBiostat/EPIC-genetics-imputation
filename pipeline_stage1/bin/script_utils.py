from __future__ import annotations

import shlex
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
    if source_build == "hg38":
        move_prefix(source_prefix, output_prefix, cwd)
        return

    triple_script = repo_root / "tools" / "triple-liftOver" / "tripleliftover_v133.pl"
    chain_gz = repo_root / "tools" / "triple-liftOver" / "library" / "chainfiles" / f"{source_build}ToHg38.over.chain.gz"
    if not triple_script.exists():
        raise FileNotFoundError(f"Missing triple-liftOver script: {triple_script}")
    if not chain_gz.exists():
        raise FileNotFoundError(f"Missing liftover chain file: {chain_gz}")
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

    run(
        plink_cmd(plink, f"--bfile {q(split_source)} --split-x b38 no-fail --make-bed --out {q(output_prefix)}"),
        cwd,
    )
    remove_prefix(split_source, cwd)

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
