#!/usr/bin/env python3
"""DevRO v5 unified CLI for structural variant discovery."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

SCRIPT_BY_MODE = {
    "dup": "VariantCaller_dup.pl",
    "inv": "VariantCaller_inv.pl",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Unified DevRO runner with publication-ready metadata, scalable genome "
            "chunking, and configurable SV calling parameters."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser("call", help="Run duplication/inversion calling.")
    call.add_argument("--mode", choices=sorted(SCRIPT_BY_MODE), required=True, help="SV mode to run.")
    call.add_argument("--regions", required=True, type=Path, help="Input region file (chr, start, chr_size).")
    call.add_argument("--config", required=True, type=Path, help="Population config (group<TAB>bam_path).")
    call.add_argument("--prefix", required=True, help="Output prefix.")
    call.add_argument("--output-dir", type=Path, default=Path("results"), help="Output directory.")
    call.add_argument("--window-size", type=int, default=1000, help="Genome scan window size.")
    call.add_argument("--read-length", type=int, default=150, help="Expected read length.")
    call.add_argument("--min-mapq", type=int, default=10, help="Minimum mapping quality.")
    call.add_argument("--mean-insert", type=float, default=400.0, help="Mean insert size.")
    call.add_argument("--insert-sigma", type=float, default=130.0, help="Insert-size standard deviation.")
    call.add_argument(
        "--chunk-lines",
        type=int,
        default=0,
        help="Split region file into chunks with this many lines (0 disables chunking).",
    )
    call.add_argument("--threads", type=int, default=1, help="Parallel chunk jobs.")

    return parser.parse_args()


def _chunk_region_file(region_path: Path, chunk_lines: int, workdir: Path) -> list[Path]:
    lines = [ln for ln in region_path.read_text().splitlines() if ln.strip()]
    if chunk_lines <= 0 or len(lines) <= chunk_lines:
        return [region_path]

    chunk_paths: list[Path] = []
    for idx in range(0, len(lines), chunk_lines):
        chunk = lines[idx : idx + chunk_lines]
        chunk_path = workdir / f"regions.chunk{idx // chunk_lines:04d}.txt"
        chunk_path.write_text("\n".join(chunk) + "\n")
        chunk_paths.append(chunk_path)
    return chunk_paths


def _run_chunk(
    perl_script: str,
    chunk_file: Path,
    prefix: str,
    config: Path,
    output_dir: Path,
    env: dict[str, str],
) -> tuple[int, list[str]]:
    cmd = ["perl", perl_script, str(chunk_file), prefix, str(config)]
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    errors = []
    if proc.returncode != 0:
        errors.append(proc.stderr.strip() or f"{perl_script} failed with exit code {proc.returncode}")
    return proc.returncode, errors


def run_call(args: argparse.Namespace) -> int:
    if shutil.which("perl") is None:
        print("ERROR: perl is required but was not found in PATH.", file=sys.stderr)
        return 2

    perl_script = SCRIPT_BY_MODE[args.mode]
    if not Path(perl_script).exists():
        print(f"ERROR: could not find {perl_script} in current directory.", file=sys.stderr)
        return 2

    if not args.regions.exists() or not args.config.exists():
        print("ERROR: --regions and --config must exist.", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_mode_dir = args.output_dir / args.mode
    output_mode_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.update(
        {
            "DEVRO_WINDOW_SIZE": str(args.window_size),
            "DEVRO_READ_LENGTH": str(args.read_length),
            "DEVRO_MIN_MAPQ": str(args.min_mapq),
            "DEVRO_MEAN_INSERT": str(args.mean_insert),
            "DEVRO_INSERT_SIGMA": str(args.insert_sigma),
            "DEVRO_OUTPUT_DIR": str(output_mode_dir),
        }
    )

    run_errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="devro-chunks-") as tmp:
        chunk_files = _chunk_region_file(args.regions, args.chunk_lines, Path(tmp))

        if len(chunk_files) == 1:
            prefix = args.prefix
            rc, errors = _run_chunk(perl_script, chunk_files[0], prefix, args.config, output_mode_dir, env)
            run_errors.extend(errors)
            if rc != 0:
                return 1
        else:
            max_workers = max(1, args.threads)
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                futures = []
                for idx, chunk in enumerate(chunk_files):
                    futures.append(
                        ex.submit(
                            _run_chunk,
                            perl_script,
                            chunk,
                            f"{args.prefix}.chunk{idx:04d}",
                            args.config,
                            output_mode_dir,
                            env,
                        )
                    )
                for fut in concurrent.futures.as_completed(futures):
                    rc, errors = fut.result()
                    run_errors.extend(errors)
                    if rc != 0:
                        return 1

    metadata = {
        "tool": "DevRO",
        "version": "5.0.0",
        "command": "call",
        "mode": args.mode,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "regions": str(args.regions),
        "config": str(args.config),
        "prefix": args.prefix,
        "output_dir": str(args.output_dir),
        "parameters": {
            "window_size": args.window_size,
            "read_length": args.read_length,
            "min_mapq": args.min_mapq,
            "mean_insert": args.mean_insert,
            "insert_sigma": args.insert_sigma,
            "chunk_lines": args.chunk_lines,
            "threads": args.threads,
        },
    }
    metadata_path = args.output_dir / f"{args.prefix}.{args.mode}.run-metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    if run_errors:
        print("\n".join(run_errors), file=sys.stderr)
        return 1

    print(f"DevRO call complete: mode={args.mode} output={output_mode_dir}")
    print(f"Run metadata written to: {metadata_path}")
    return 0


def main() -> int:
    args = parse_args()
    if args.command == "call":
        return run_call(args)
    print("Unsupported command", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
