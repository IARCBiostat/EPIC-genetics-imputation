#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "pipeline_stage1" / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from stage1_reports import TASK_SEQUENCE, run_task  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create stage-1 figures, tables, flags, and report assets.")
    parser.add_argument("--analysis-root", default=str(REPO_ROOT / "analysis"))
    parser.add_argument("--studies", default="all")
    parser.add_argument("--reference-dir", default=str(REPO_ROOT / "data" / "reference" / "1000G"))
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()
    reference_dir = Path(args.reference_dir).resolve()
    for task_id in TASK_SEQUENCE:
        run_task(task_id, analysis_root, args.studies, reference_dir, args.force)


if __name__ == "__main__":
    main()
