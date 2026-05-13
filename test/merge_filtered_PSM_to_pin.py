#!/usr/bin/env python3
"""Merge Percolator filtered PSMs with Sipros .pin feature columns.

This is a standalone version of script33.assembly.intergrate_filtered_psms_with_feature().
It reads target/decoy Percolator PSM tables, filters them at 1% q-value, and joins
the selected PSMs back to the Sipros .pin file on PSMId == SpecId.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_SEARCH_DIR = Path("data/search")
DEFAULT_SAMPLE = "Pan_062822_X1iso5"
PERCOLATOR_COLUMNS = ["PSMId", "score", "q-value", "posterior_error_prob"]


def default_output_path(pin_path: Path) -> Path:
    """Return the script33-style *_filtered_psms.tsv output path."""
    if pin_path.suffix == ".pin":
        return pin_path.with_name(f"{pin_path.stem}_filtered_psms.tsv")
    return pin_path.with_suffix(".filtered_psms.tsv")


def sample_paths(search_dir: Path, sample: str) -> tuple[Path, Path, Path]:
    """Return target, decoy, and pin paths for a sample name."""
    return (
        search_dir / f"{sample}_target_psms.tsv",
        search_dir / f"{sample}_decoy_psms.tsv",
        search_dir / f"{sample}.pin",
    )


def read_percolator_psms(path: Path, q_value_threshold: float) -> pd.DataFrame:
    """Read and q-value filter a Percolator PSM TSV."""
    psm = pd.read_csv(path, sep="\t")
    missing = [column for column in PERCOLATOR_COLUMNS if column not in psm.columns]
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")
    return psm.loc[psm["q-value"] <= q_value_threshold, PERCOLATOR_COLUMNS]


def merge_filtered_psms(
    target_psms_path: Path,
    decoy_psms_path: Path,
    pin_path: Path,
    output_path: Path,
    q_value_threshold: float = 0.01,
) -> pd.DataFrame:
    """Merge filtered target/decoy Percolator PSMs with .pin feature data."""
    filtered_target = read_percolator_psms(target_psms_path, q_value_threshold)
    filtered_decoy = read_percolator_psms(decoy_psms_path, q_value_threshold)
    psm = pd.concat([filtered_target, filtered_decoy], ignore_index=True)

    pin = pd.read_csv(pin_path, sep="\t")
    if "SpecId" not in pin.columns:
        raise ValueError(f"{pin_path} is missing required column: SpecId")

    merged = pd.merge(psm, pin, left_on="PSMId", right_on="SpecId", how="left")
    missing_pin_matches = int(merged["SpecId"].isna().sum())
    if missing_pin_matches:
        raise ValueError(
            f"{missing_pin_matches} filtered PSM(s) from Percolator did not match "
            f"any SpecId in {pin_path}"
        )

    merged = merged.drop(columns=["SpecId"])
    if "ScanNr" in merged.columns:
        merged["ScanNr"] = merged["ScanNr"].astype(int)
        merged = merged.sort_values(by="ScanNr")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter target/decoy Percolator PSMs by q-value and merge them with "
            "Sipros .pin feature columns."
        )
    )
    parser.add_argument(
        "sample",
        nargs="?",
        default=DEFAULT_SAMPLE,
        help=(
            "Sample basename used to derive input paths as "
            "<search-dir>/<sample>_target_psms.tsv, "
            "<search-dir>/<sample>_decoy_psms.tsv, and <search-dir>/<sample>.pin. "
            f"Defaults to {DEFAULT_SAMPLE}."
        ),
    )
    parser.add_argument(
        "--search-dir",
        type=Path,
        default=DEFAULT_SEARCH_DIR,
        help=f"Directory containing sample Percolator and .pin files. Defaults to {DEFAULT_SEARCH_DIR}.",
    )
    parser.add_argument(
        "--target-psms",
        type=Path,
        default=None,
        help="Percolator target PSM TSV. Overrides the path derived from sample.",
    )
    parser.add_argument(
        "--decoy-psms",
        type=Path,
        default=None,
        help="Percolator decoy PSM TSV. Overrides the path derived from sample.",
    )
    parser.add_argument(
        "--pin",
        type=Path,
        default=None,
        help="Sipros .pin file containing feature columns. Overrides the path derived from sample.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output TSV path. Defaults to <pin_stem>_filtered_psms.tsv.",
    )
    parser.add_argument(
        "--q-value",
        type=float,
        default=0.01,
        help="q-value threshold used to keep Percolator PSMs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default_target_psms, default_decoy_psms, default_pin = sample_paths(args.search_dir, args.sample)
    target_psms = args.target_psms or default_target_psms
    decoy_psms = args.decoy_psms or default_decoy_psms
    pin = args.pin or default_pin
    output = args.output or default_output_path(pin)
    merged = merge_filtered_psms(
        target_psms_path=target_psms,
        decoy_psms_path=decoy_psms,
        pin_path=pin,
        output_path=output,
        q_value_threshold=args.q_value,
    )
    print(f"Wrote {len(merged)} filtered PSMs to {output}")


if __name__ == "__main__":
    main()
