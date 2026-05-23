#!/usr/bin/env python3
"""Run pct1-100 FT2 searches, Percolator filtering, and PSM/PIN merging.

The default benchmark processes the main WinnowNet pct folders:
../WinnowNet/data/{pct1,pct2,pct5,pct25,pct50,pct99}.
Each FT2 file is searched against data/spectra, filtered with Percolator, and
merged into data/search/<sample>_filtered_psms.tsv.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from merge_filtered_PSM_to_pin import (
    default_output_path,
    merge_filtered_psms,
    sample_paths,
)


DEFAULT_PCT_DIRS = ("pct1", "pct2", "pct5", "pct25", "pct50", "pct99")
DEFAULT_WINNOWNET_DATA = Path("../WinnowNet/data")
DEFAULT_SPECTRA_DIR = Path("data/spectra")
DEFAULT_SEARCH_DIR = Path("data/search")
DEFAULT_CONFIG = Path("configTemplates/SIP.cfg")
DEFAULT_SEARCH_TOOL = Path("tools/sipros_search_spectra")
DEFAULT_PERCOLATOR = Path("tools/percolator")
DEFAULT_THREADS_PER_FILE = 8
DEFAULT_TOTAL_THREADS = 24
DEFAULT_Q_VALUE = 0.01


@dataclass(frozen=True)
class BenchmarkConfig:
    repo_root: Path
    winnownet_data: Path
    spectra_dir: Path
    search_dir: Path
    config: Path
    search_tool: Path
    percolator: Path
    threads_per_file: int
    q_value: float
    dry_run: bool


@dataclass(frozen=True)
class Ft2Job:
    sample: str
    ft2_path: Path
    ft1_path: Path


@dataclass(frozen=True)
class JobResult:
    sample: str
    pin_path: Path
    target_psms_path: Path
    decoy_psms_path: Path
    filtered_psms_path: Path
    filtered_psms: int | None


def existing_ft1_path(ft2_path: Path) -> Path | None:
    """Return the matching FT1 path next to an FT2 file, if present."""
    candidates = [
        ft2_path.with_suffix(".FT1"),
        ft2_path.with_suffix(".ft1"),
        ft2_path.with_suffix(".Ft1"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def discover_ft2_jobs(
    winnownet_data: Path,
    pct_dirs: Sequence[str],
    samples: set[str] | None,
) -> list[Ft2Job]:
    """Find FT2 files in the selected pct directories."""
    jobs: list[Ft2Job] = []
    for pct_dir in pct_dirs:
        sample_dir = winnownet_data / pct_dir
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"Missing pct directory: {sample_dir}")

        ft2_paths = sorted(
            path
            for pattern in ("*.FT2", "*.ft2")
            for path in sample_dir.glob(pattern)
            if path.is_file()
        )
        for ft2_path in ft2_paths:
            sample = ft2_path.stem
            if samples is not None and sample not in samples:
                continue

            ft1_path = existing_ft1_path(ft2_path)
            if ft1_path is None:
                raise FileNotFoundError(f"No matching FT1 file found for {ft2_path}")
            jobs.append(Ft2Job(sample=sample, ft2_path=ft2_path, ft1_path=ft1_path))

    if not jobs:
        selected = ", ".join(pct_dirs)
        if samples:
            selected += f"; samples={', '.join(sorted(samples))}"
        raise FileNotFoundError(f"No FT2 files found for {selected}")

    seen: dict[str, Path] = {}
    for job in jobs:
        previous = seen.get(job.sample)
        if previous is not None:
            raise ValueError(
                f"Duplicate FT2 basename {job.sample!r}: {previous} and {job.ft2_path}. "
                "Use a narrower --pct-dirs/--samples selection or unique output naming."
            )
        seen[job.sample] = job.ft2_path

    return jobs


def validate_inputs(config: BenchmarkConfig) -> None:
    """Check required files and directories before launching workers."""
    required_dirs = {
        "WinnowNet data directory": config.winnownet_data,
        "spectra directory": config.spectra_dir,
    }
    for label, path in required_dirs.items():
        if not path.is_dir():
            raise FileNotFoundError(f"Missing {label}: {path}")

    required_files = {
        "Sipros config": config.config,
        "sipros_search_spectra": config.search_tool,
        "Percolator": config.percolator,
    }
    for label, path in required_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")

    if config.threads_per_file < 1:
        raise ValueError("--threads-per-file must be at least 1")
    if config.q_value < 0:
        raise ValueError("--q-value must be non-negative")


def max_workers(total_threads: int, threads_per_file: int) -> int:
    if total_threads < 1:
        raise ValueError("--total-threads must be at least 1")
    if threads_per_file < 1:
        raise ValueError("--threads-per-file must be at least 1")
    return max(1, total_threads // threads_per_file)


def format_command(command: Sequence[object]) -> str:
    """Return a readable shell-like command for logging."""
    return " ".join(str(part) for part in command)


def search_command(job: Ft2Job, config: BenchmarkConfig) -> list[str]:
    return [
        str(config.search_tool),
        "-f",
        str(job.ft2_path),
        "-c",
        str(config.config),
        "-h5",
        str(config.spectra_dir),
        "-o",
        str(config.search_dir),
        "--rt-tolerance",
        "5",
        "--tolerance",
        "10",
        "--tolerance-unit",
        "ppm",
        "-t",
        str(config.threads_per_file),
    ]


def percolator_command(
    pin_path: Path,
    target_psms_path: Path,
    decoy_psms_path: Path,
    config: BenchmarkConfig,
) -> list[str]:
    return [
        str(config.percolator),
        "--only-psms",
        "--no-terminate",
        "--num-threads",
        str(config.threads_per_file),
        "--results-psms",
        str(target_psms_path),
        "--decoy-results-psms",
        str(decoy_psms_path),
        str(pin_path),
    ]


def run_command(command: Sequence[str], config: BenchmarkConfig, sample: str, stage: str) -> None:
    """Run a subprocess and raise a clear error if it fails."""
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(config.threads_per_file)

    try:
        completed = subprocess.run(
            command,
            cwd=config.repo_root,
            env=env,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        message = [
            f"{stage} failed for {sample} with exit code {exc.returncode}",
            f"Command: {format_command(command)}",
        ]
        if exc.stdout:
            message.append(f"stdout:\n{exc.stdout.rstrip()}")
        if exc.stderr:
            message.append(f"stderr:\n{exc.stderr.rstrip()}")
        raise RuntimeError("\n".join(message)) from exc

    if completed.stdout:
        print(completed.stdout.rstrip())
    if completed.stderr:
        print(completed.stderr.rstrip())


def require_output(path: Path, sample: str, stage: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{stage} did not produce expected output for {sample}: {path}")


def process_job(job: Ft2Job, config: BenchmarkConfig) -> JobResult:
    """Search one FT2 file, run Percolator, and merge filtered PSMs."""
    target_psms_path, decoy_psms_path, pin_path = sample_paths(config.search_dir, job.sample)
    filtered_psms_path = default_output_path(pin_path)

    search_cmd = search_command(job, config)
    percolator_cmd = percolator_command(pin_path, target_psms_path, decoy_psms_path, config)

    if config.dry_run:
        print(f"[dry-run] {job.sample}")
        print(f"  FT2: {job.ft2_path}")
        print(f"  FT1: {job.ft1_path}")
        print(f"  search: {format_command(search_cmd)}")
        print(f"  percolator: {format_command(percolator_cmd)}")
        print(
            "  merge: "
            f"{target_psms_path} + {decoy_psms_path} + {pin_path} -> {filtered_psms_path}"
        )
        return JobResult(
            sample=job.sample,
            pin_path=pin_path,
            target_psms_path=target_psms_path,
            decoy_psms_path=decoy_psms_path,
            filtered_psms_path=filtered_psms_path,
            filtered_psms=None,
        )

    config.search_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{job.sample}] search")
    run_command(search_cmd, config, job.sample, "sipros_search_spectra")
    require_output(pin_path, job.sample, "sipros_search_spectra")

    print(f"[{job.sample}] percolator")
    run_command(percolator_cmd, config, job.sample, "Percolator")
    require_output(target_psms_path, job.sample, "Percolator")
    require_output(decoy_psms_path, job.sample, "Percolator")

    print(f"[{job.sample}] merge")
    merged = merge_filtered_psms(
        target_psms_path=target_psms_path,
        decoy_psms_path=decoy_psms_path,
        pin_path=pin_path,
        output_path=filtered_psms_path,
        q_value_threshold=config.q_value,
    )
    require_output(filtered_psms_path, job.sample, "merge")

    return JobResult(
        sample=job.sample,
        pin_path=pin_path,
        target_psms_path=target_psms_path,
        decoy_psms_path=decoy_psms_path,
        filtered_psms_path=filtered_psms_path,
        filtered_psms=len(merged),
    )


def process_jobs(jobs: Sequence[Ft2Job], config: BenchmarkConfig, workers: int) -> list[JobResult]:
    """Run FT2 jobs with fail-fast exception handling."""
    if config.dry_run or workers == 1:
        return [process_job(job, config) for job in jobs]

    results: list[JobResult] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_job = {executor.submit(process_job, job, config): job for job in jobs}
        for future in concurrent.futures.as_completed(future_to_job):
            job = future_to_job[future]
            try:
                result = future.result()
            except Exception:
                for pending in future_to_job:
                    pending.cancel()
                executor.shutdown(wait=False, cancel_futures=True)
                raise
            results.append(result)
            print(f"[{job.sample}] wrote {result.filtered_psms_path}")
    return sorted(results, key=lambda result: result.sample)


def comma_split(values: Iterable[str]) -> list[str]:
    items: list[str] = []
    for value in values:
        items.extend(part.strip() for part in value.split(",") if part.strip())
    return items


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search WinnowNet pct FT2 files against SIP spectra, run Percolator, "
            "and merge filtered PSMs back to PIN features."
        )
    )
    parser.add_argument(
        "--winnownet-data",
        type=Path,
        default=DEFAULT_WINNOWNET_DATA,
        help=f"WinnowNet data root. Defaults to {DEFAULT_WINNOWNET_DATA}.",
    )
    parser.add_argument(
        "--pct-dirs",
        nargs="+",
        default=list(DEFAULT_PCT_DIRS),
        help=(
            "pct directories under --winnownet-data. Accepts spaces or comma-separated "
            f"values. Defaults to {', '.join(DEFAULT_PCT_DIRS)}."
        ),
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Optional FT2 basenames to process. Accepts spaces or comma-separated values.",
    )
    parser.add_argument(
        "--spectra-dir",
        type=Path,
        default=DEFAULT_SPECTRA_DIR,
        help=f"SIP spectra HDF5 directory. Defaults to {DEFAULT_SPECTRA_DIR}.",
    )
    parser.add_argument(
        "--search-dir",
        type=Path,
        default=DEFAULT_SEARCH_DIR,
        help=(
            "Directory for .pin, Percolator TSV, and merged TSV outputs. "
            f"Defaults to {DEFAULT_SEARCH_DIR}."
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"Sipros config file. Defaults to {DEFAULT_CONFIG}.",
    )
    parser.add_argument(
        "--search-tool",
        type=Path,
        default=DEFAULT_SEARCH_TOOL,
        help=f"sipros_search_spectra executable. Defaults to {DEFAULT_SEARCH_TOOL}.",
    )
    parser.add_argument(
        "--percolator",
        type=Path,
        default=DEFAULT_PERCOLATOR,
        help=f"Percolator executable. Defaults to {DEFAULT_PERCOLATOR}.",
    )
    parser.add_argument(
        "--threads-per-file",
        type=int,
        default=DEFAULT_THREADS_PER_FILE,
        help=(
            "Threads used by each FT2 search and Percolator run. "
            f"Defaults to {DEFAULT_THREADS_PER_FILE}."
        ),
    )
    parser.add_argument(
        "--total-threads",
        type=int,
        default=DEFAULT_TOTAL_THREADS,
        help=f"Total thread budget. Defaults to {DEFAULT_TOTAL_THREADS}.",
    )
    parser.add_argument(
        "--q-value",
        type=float,
        default=DEFAULT_Q_VALUE,
        help=f"q-value threshold for merged filtered PSMs. Defaults to {DEFAULT_Q_VALUE}.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print discovered jobs and commands without running search, Percolator, or merge.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path.cwd()
    pct_dirs = comma_split(args.pct_dirs)
    samples = set(comma_split(args.samples)) if args.samples else None

    config = BenchmarkConfig(
        repo_root=repo_root,
        winnownet_data=args.winnownet_data,
        spectra_dir=args.spectra_dir,
        search_dir=args.search_dir,
        config=args.config,
        search_tool=args.search_tool,
        percolator=args.percolator,
        threads_per_file=args.threads_per_file,
        q_value=args.q_value,
        dry_run=args.dry_run,
    )

    validate_inputs(config)
    jobs = discover_ft2_jobs(config.winnownet_data, pct_dirs, samples)
    workers = max_workers(args.total_threads, args.threads_per_file)

    print(
        f"Processing {len(jobs)} FT2 file(s) with {workers} worker(s), "
        f"{args.threads_per_file} thread(s) per file."
    )
    results = process_jobs(jobs, config, workers)

    if args.dry_run:
        print(f"Dry run listed {len(results)} job(s).")
        return

    print("Completed filtered PSM outputs:")
    for result in results:
        print(f"  {result.sample}: {result.filtered_psms_path} ({result.filtered_psms} rows)")


if __name__ == "__main__":
    main()
