#!/usr/bin/env python3
"""Concatenate Pfam HMM profiles and build a pressed HMMER database."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import traceback
from pathlib import Path


def resolve_path(pfam_path: str, pipeline_root: Path) -> Path:
    """Resolve a Pfam HMM file path relative to common search locations.

    Args:
        pfam_path: User-supplied path to a Pfam HMM file (absolute or relative).
        pipeline_root: Root directory of the finder pipeline.

    Returns:
        Resolved absolute path to an existing Pfam HMM file.

    Raises:
        FileNotFoundError: If the file cannot be found in any searched location.
    """
    candidate = Path(pfam_path)
    if candidate.is_file():
        return candidate

    root_candidate = pipeline_root / pfam_path
    if root_candidate.is_file():
        return root_candidate

    cwd_candidate = Path.cwd() / pfam_path
    if cwd_candidate.is_file():
        return cwd_candidate

    raise FileNotFoundError(f"Pfam HMM file not found: {pfam_path}")


def main() -> None:
    """Combine Pfam HMM files and run ``hmmpress`` on the merged database."""
    parser = argparse.ArgumentParser(
        description="Combine Pfam HMM profiles and run hmmpress."
    )
    parser.add_argument("--output", required=True, help="Output combined HMM path.")
    parser.add_argument("--log", required=True, help="Log file path.")
    parser.add_argument(
        "--pfam",
        dest="pfams",
        action="append",
        required=True,
        help="Pfam HMM file; repeat for multiple profiles.",
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    log_path = Path(args.log)
    pipeline_root = Path(__file__).resolve().parents[1]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with log_path.open("w") as log:
            log.write("Starting build_combined_hmm\n")
            log.write(f"Output: {output_path}\n")
            log.write(f"Pfam files: {args.pfams}\n")
            log.write(f"Pipeline root: {pipeline_root}\n")

            resolved_pfams = []
            for pfam in args.pfams:
                resolved = resolve_path(pfam, pipeline_root)
                resolved_pfams.append(resolved)
                log.write(f"Resolved {pfam} -> {resolved}\n")

            log.write(f"Combining {len(resolved_pfams)} Pfam files...\n")
            with output_path.open("w") as out_f:
                for pfam_path in resolved_pfams:
                    log.write(f"Reading {pfam_path}\n")
                    with pfam_path.open() as pfam_f:
                        out_f.write(pfam_f.read())

            log.write(f"Combined {len(resolved_pfams)} Pfam files into {output_path}\n")
            log.write(f"Running hmmpress on {output_path}...\n")
            for suffix in (".h3f", ".h3p", ".h3i", ".h3m"):
                extra = Path(f"{output_path}{suffix}")
                if extra.exists():
                    extra.unlink()
            result = subprocess.run(
                ["hmmpress", str(output_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
            log.write(result.stdout or "")
            if result.stderr:
                log.write(result.stderr)
            log.write("hmmpress completed successfully\n")

    except Exception as exc:
        error_msg = f"Error in build_combined_hmm: {type(exc).__name__}: {exc}\n"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a") as log:
            log.write(error_msg)
            log.write(traceback.format_exc())
        print(error_msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
