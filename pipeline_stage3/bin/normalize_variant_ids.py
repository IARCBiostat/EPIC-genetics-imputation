#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prefer rsIDs and fall back to chr:pos:ref:alt IDs.")
    parser.add_argument("--input", required=True, help="Input VCF or VCF.GZ file.")
    parser.add_argument("--output", required=True, help="Output VCF or VCF.GZ file.")
    parser.add_argument("--mapping-output", required=True, help="Output TSV/TSV.GZ mapping file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    mapping_path = Path(args.mapping_output)

    seen_ids: set[str] = set()

    with open_text(input_path, "r") as src, open_text(output_path, "w") as dst, open_text(mapping_path, "w") as mapping:
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
