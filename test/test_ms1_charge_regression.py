#!/usr/bin/env python3
"""Regression checks for Raxport HDF5-only scan input.

The synthetic files use the schema v5 groups/datasets read by Sipros. The test
keeps the data tiny so it can run after a local build without relying on raw
vendor fixtures.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "data" / "tmp" / "raxport_hdf5_workflow_test"
SIPROS = ROOT / "bin" / "sipros"
CFG = ROOT / "configTemplates" / "Regular.cfg"
FASTA = ROOT / "data" / "EcoliWithCrapNodup.fasta"


def write_dataset(group: h5py.Group, name: str, values, dtype) -> None:
    group.create_dataset(name, data=np.asarray(values, dtype=dtype))


def create_raxport_hdf5(path: Path, ms_orders: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    peak_mz: list[float] = []
    peak_intensity: list[float] = []
    peak_charge: list[int] = []
    peak_start: list[int] = []
    peak_count: list[int] = []
    reaction_start: list[int] = []
    reaction_count: list[int] = []
    parent_scan_number: list[int] = []
    precursor_mass: list[float] = []
    reaction_charge: list[int] = []
    candidate_start: list[int] = []
    candidate_count: list[int] = []
    candidate_charge: list[int] = []
    candidate_mz: list[float] = []
    candidate_intensity: list[float] = []

    for idx, order in enumerate(ms_orders):
        scan_number = 1001 + idx
        if order == 1:
            mz = [499.8, 500.3, 501.2]
            intensity = [1000.0, 2000.0, 500.0]
        else:
            mz = [110.0, 175.0, 240.0, 310.0]
            intensity = [50.0, 75.0, 125.0, 25.0]
        peak_start.append(len(peak_mz))
        peak_count.append(len(mz))
        peak_mz.extend(mz)
        peak_intensity.extend(intensity)
        peak_charge.extend([0] * len(mz))

        if order >= 2:
            reaction_start.append(len(precursor_mass))
            reaction_count.append(1)
            parent_scan_number.append(1001)
            precursor_mass.append(500.3)
            reaction_charge.append(2)
            candidate_start.append(len(candidate_mz))
            candidate_count.append(1)
            candidate_charge.append(2)
            candidate_mz.append(500.3)
            candidate_intensity.append(2000.0)
        else:
            reaction_start.append(-1)
            reaction_count.append(0)
            parent_scan_number.append(0)

    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = np.int32(5)
        handle.attrs["source_raw_file"] = b"synthetic.raw"
        handle.attrs["instrument_model"] = b"synthetic"
        handle.attrs["raxport_version"] = b"test"

        scans = handle.create_group("scans")
        write_dataset(scans, "scan_number", [1001 + i for i in range(len(ms_orders))], np.int32)
        write_dataset(scans, "ms_order", ms_orders, np.int32)
        write_dataset(scans, "retention_time", [1.0 + i * 0.1 for i in range(len(ms_orders))], np.float64)
        write_dataset(scans, "tic", [sum(peak_intensity)] * len(ms_orders), np.float64)
        write_dataset(scans, "scan_filter_id", [0] * len(ms_orders), np.int32)
        write_dataset(scans, "activation_id", [0] * len(ms_orders), np.int32)
        write_dataset(scans, "parent_scan_number", parent_scan_number, np.int32)
        write_dataset(scans, "reaction_start", reaction_start, np.int64)
        write_dataset(scans, "reaction_count", reaction_count, np.int32)
        write_dataset(scans, "peak_start", peak_start, np.int64)
        write_dataset(scans, "peak_count", peak_count, np.int32)

        peaks = handle.create_group("peaks")
        write_dataset(peaks, "mz", peak_mz, np.float64)
        write_dataset(peaks, "intensity", peak_intensity, np.float64)
        write_dataset(peaks, "resolution", [0.0] * len(peak_mz), np.float64)
        write_dataset(peaks, "baseline", [0.0] * len(peak_mz), np.float64)
        write_dataset(peaks, "noise", [0.0] * len(peak_mz), np.float64)
        write_dataset(peaks, "charge", peak_charge, np.int32)
        write_dataset(peaks, "mobility_trace_start", [-1] * len(peak_mz), np.int64)
        write_dataset(peaks, "mobility_trace_count", [0] * len(peak_mz), np.int32)

        reactions = handle.create_group("reactions")
        n_reactions = len(precursor_mass)
        write_dataset(reactions, "precursor_mass", precursor_mass, np.float64)
        write_dataset(reactions, "isolation_width", [1.0] * n_reactions, np.float64)
        write_dataset(reactions, "charge_state", reaction_charge, np.int32)
        write_dataset(reactions, "collision_energy", [0.0] * n_reactions, np.float64)
        write_dataset(reactions, "collision_energy_valid", [0] * n_reactions, np.int32)
        write_dataset(reactions, "activation_type_id", [0] * n_reactions, np.int32)
        write_dataset(reactions, "multiple_activation", [0] * n_reactions, np.int32)
        write_dataset(reactions, "precursor_range_valid", [0] * n_reactions, np.int32)
        write_dataset(reactions, "first_precursor_mass", precursor_mass, np.float64)
        write_dataset(reactions, "last_precursor_mass", precursor_mass, np.float64)
        write_dataset(reactions, "isolation_width_offset", [0.0] * n_reactions, np.float64)
        write_dataset(reactions, "one_over_k0_begin", [0.0] * n_reactions, np.float64)
        write_dataset(reactions, "one_over_k0_end", [0.0] * n_reactions, np.float64)
        write_dataset(reactions, "candidate_start", candidate_start, np.int64)
        write_dataset(reactions, "candidate_count", candidate_count, np.int32)

        candidates = handle.create_group("precursor_candidates")
        write_dataset(candidates, "charge", candidate_charge, np.int32)
        write_dataset(candidates, "mz", candidate_mz, np.float64)
        write_dataset(candidates, "intensity", candidate_intensity, np.float64)
        write_dataset(candidates, "one_over_k0", [0.0] * len(candidate_mz), np.float64)

        traces = handle.create_group("peak_mobility_traces")
        write_dataset(traces, "one_over_k0_index", [], np.int32)
        write_dataset(traces, "intensity", [], np.float32)

        strings = handle.create_group("string_tables")
        strings.create_dataset("scan_filter", data=np.asarray([b"synthetic"], dtype="S64"))
        strings.create_dataset("activation", data=np.asarray([b"none"], dtype="S64"))
        strings.create_dataset("reaction_activation_type", data=np.asarray([b"none"], dtype="S64"))


def run_sipros(scan_file: Path, out_dir: Path) -> subprocess.CompletedProcess[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        [
            str(SIPROS),
            "search-fasta",
            "-f", str(scan_file),
            "-c", str(CFG),
            "-fasta", str(FASTA),
            "-o", str(out_dir),
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> None:
    if not SIPROS.exists():
        raise SystemExit(f"Build sipros first: {SIPROS}")
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)

    bad_input = TMP / "bad.scan"
    bad_input.write_text("S\t1\t500.0\n", encoding="utf-8")
    bad = run_sipros(bad_input, TMP / "bad_out")
    assert bad.returncode != 0
    assert "Raxport HDF5 scan input required" in bad.stdout

    mixed = TMP / "mixed_ms2_ms3.h5"
    create_raxport_hdf5(mixed, [1, 2, 3])
    result = run_sipros(mixed, TMP / "mixed_out")
    assert result.returncode == 0, result.stdout
    assert "Preprocessing scans: 1" in result.stdout, result.stdout

    ms3_only = TMP / "ms3_only.h5"
    create_raxport_hdf5(ms3_only, [1, 3])
    result = run_sipros(ms3_only, TMP / "ms3_out")
    assert result.returncode != 0, result.stdout
    assert "No ms_order == 2 scans found" in result.stdout, result.stdout

    print("ok: Raxport HDF5-only input rejects non-HDF5 input, loads ms_order 2, and ignores ms_order > 2")


if __name__ == "__main__":
    main()
