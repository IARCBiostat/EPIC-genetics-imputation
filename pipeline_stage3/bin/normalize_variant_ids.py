#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


@contextmanager
def open_text(path: str, mode: str) -> Iterator[TextIO]:
    if path == "-":
        if "r" in mode:
            yield sys.stdin
        else:
            yield sys.stdout
        return

    file_path = Path(path)
    if file_path.suffix == ".gz":
        with gzip.open(file_path, mode + "t") as handle:
            yield handle
        return

    with file_path.open(mode) as handle:
        yield handle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prefer rsIDs and fall back to chr:pos:ref:alt IDs.")
    parser.add_argument("--input", required=True, help="Input VCF or VCF.GZ file.")
    parser.add_argument("--output", required=True, help="Output VCF or VCF.GZ file.")
    parser.add_argument("--mapping-output", required=True, help="Output TSV/TSV.GZ mapping file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    seen_ids: set[str] = set()

    with open_text(args.input, "r") as src, open_text(args.output, "w") as dst, open_text(args.mapping_output, "w") as mapping:
        mapping.write("CHROM\tPOS\tREF\tALT\tORIGINAL_ID\tFINAL_ID\tID_SOURCE\n")

        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom, pos, original_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            fallback_id = f"{chrom}:{pos}:{ref}:{alt}"
            final_id = fallback_id
            source = "fallback"

            if original_id and original_id != "." and original_id.lower().startswith("rs"):
                if original_id not in seen_ids:
                    final_id = original_id
                    source = "rsid"
                    seen_ids.add(original_id)
                else:
                    source = "duplicate_rsid_fallback"

            fields[2] = final_id
            dst.write("\t".join(fields) + "\n")
            mapping.write(
                "\t".join(
                    [
                        chrom,
                        pos,
                        ref,
                        alt,
                        original_id,
                        final_id,
                        source,
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
