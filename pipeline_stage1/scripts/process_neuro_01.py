#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))

from script_utils import (
    configure_plink,
    apply_epic_sample_metadata,
    ensure_dirs,
    liftover_to_hg38,
    plink_cmd,
    q,
    remove_prefix,
    run,
    write_summary,
)


PREFIX = 'Neuro_01'
STUDY_ID = 'Neuro_01'
RAW_REL = 'genetics/Neuro_01/Data_Received/Genetic data before QC/GSA_epic'
RAW_INPUT_MODE = 'ped_map'
MANIFEST_REL = None
BUILD = '37'
ID_LINK_REL = None


def main() -> None:
    parser = argparse.ArgumentParser(description=f"Stage 1 processing for {STUDY_ID}")
    parser.add_argument('--data-root', required=True, help='Root directory containing the current EPIC genetics data layout.')
    parser.add_argument('--work-root', default=str(Path(__file__).resolve().parents[1] / 'work'), help='Output root for stage-1 processing results.')
    parser.add_argument('--outdir', default=str(Path(__file__).resolve().parents[2] / 'analysis'), help='Analysis root for final stage-1 outputs.')
    parser.add_argument('--plink', default='plink', help='PLINK executable to use.')
    parser.add_argument('--python2', default='python2.7', help='Unused for this study; kept for wrapper compatibility.')
    parser.add_argument(
        '--build',
        default=os.environ.get('NEURO_01_BUILD', BUILD),
        choices=['36', '37', '38'],
        help='Source genome build for the delivered Neuro_01 PED/MAP dataset.',
    )
    args = parser.parse_args()
    args.plink = configure_plink(args.plink)

    data_root = Path(args.data_root).resolve()
    work_root = Path(args.work_root).resolve()
    pipeline_root = Path(__file__).resolve().parents[1]
    study_dir = work_root / STUDY_ID
    trace = study_dir / 'Trace'
    lift_dir = trace / 'LiftOver'
    stage1_dir = Path(args.outdir).resolve() / STUDY_ID / 'stage1'
    summary = stage1_dir / 'summary.txt'

    ensure_dirs([study_dir, trace, lift_dir, stage1_dir])

    base_prefix = trace / PREFIX
    metadata_prefix = trace / f'{PREFIX}_metadata'
    ready_prefix = trace / f'{PREFIX}_ready'
    final_prefix = stage1_dir / STUDY_ID

    run(
        plink_cmd(args.plink, f"--file {q(data_root / RAW_REL)} --make-bed --out {q(base_prefix)}"),
        study_dir,
    )

    write_summary(
        summary,
        (
            "\n*************************** Neuro_01 Stage 1 *******************************\n"
            "Using the delivered before-QC PED/MAP dataset as the stage-1 input.\n"
            "Converted the PED/MAP delivery to a PLINK binary dataset before stage-1 processing.\n"
            "No study-specific chip manifest or EPIC linkage file is available in the current synced layout,\n"
            "so manifest-driven normalization and EPIC ID remapping are skipped for this study.\n"
            f"Assumed source build: {args.build}\n\n"
        ),
        'w',
    )

    apply_epic_sample_metadata(base_prefix, metadata_prefix, data_root, STUDY_ID, summary, strict=False)
    run(
        plink_cmd(args.plink, f"--bfile {q(metadata_prefix)} --indiv-sort n --make-bed --out {q(ready_prefix)}"),
        study_dir,
    )
    remove_prefix(metadata_prefix, study_dir)
    remove_prefix(base_prefix, study_dir)

    write_summary(
        summary,
        (
            "\n*************************** Pre-Liftover Preparation *******************************\n"
            "Assigned sex and phenotype from the shared EPIC study case-status file.\n"
            "Sorted samples by individual ID before liftover.\n\n"
        ),
        'a',
    )

    liftover_to_hg38(
        args.plink,
        pipeline_root.parent,
        ready_prefix,
        final_prefix,
        args.build,
        lift_dir,
        summary,
        study_dir,
    )
    remove_prefix(ready_prefix, study_dir)


if __name__ == "__main__":
    main()
