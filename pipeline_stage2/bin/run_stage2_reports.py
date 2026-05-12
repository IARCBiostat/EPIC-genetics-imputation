#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "pipeline_stage2" / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from stage2_artifacts import run_reports  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Stage 2 summary and per-study report artifacts.")
    parser.add_argument("--analysis-root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--stage1-root", default=None)
    parser.add_argument("--summary-output", default=None)
    parser.add_argument("--studies", default="all")
    parser.add_argument("--chromosomes", default="all")
    parser.add_argument("--reference-dir", default=str(REPO_ROOT / "data" / "reference" / "1000G"))
    parser.add_argument("--min-r2-threshold", default="0.3")
    parser.add_argument("--high-quality-threshold", default="0.8")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    stage1_root = Path(args.stage1_root).resolve() if args.stage1_root else analysis_root
    summary_output = Path(args.summary_output).resolve() if args.summary_output else analysis_root / "stage2-summary.md"
    reference_dir = Path(args.reference_dir).resolve()

    run_reports(analysis_root, args.studies, reference_dir, args.force, args.chromosomes)


if __name__ == "__main__":
    main()
