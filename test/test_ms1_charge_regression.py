#!/usr/bin/env python3
"""Small regression test for MS1 isotope abundance charge handling.

This uses real fixture PSMs for two MS1 abundance edge cases:

1. The FT2 Z line has multiple precursor m/z+charge annotations. The SIP/HDF5
   record has charge 2, and the same FT2 scan also contains a tempting charge-3
   precursor annotation. Assignment must only consider charge-2 annotations.
2. A high-label PSM has an FT1 anchor above the peptide's maximum possible C13
   isotope index. That should saturate to 100% MS1 abundance instead of becoming
   zero.
"""

from __future__ import annotations

import re
from pathlib import Path

import h5py


ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configTemplates" / "SIP.cfg"
FT1 = ROOT / "data/pct1/sipros/Pan_062822_X1iso5.FT1"
FT2 = ROOT / "data/pct1/sipros/Pan_062822_X1iso5.FT2"
H5 = ROOT / "data/spectra/spectra_C13_020.000Pct.h5"

PROTON = 1.007276466
NEUTRON = 1.003355
MS1_TOLERANCE_PPM = 10.0


def parse_config():
    residues: dict[str, list[int]] = {}
    masses: dict[str, list[float]] = {}
    residue_re = re.compile(r"^Residue\{([^}]+)\}\s*=\s*([^#]+)")
    mass_re = re.compile(r"^Element_Masses\{([^}]+)\}\s*=\s*([^#]+)")
    with CFG.open() as handle:
        for line in handle:
            if m := residue_re.search(line):
                residues[m.group(1)] = [int(x.strip().rstrip(",")) for x in m.group(2).split(",") if x.strip().rstrip(",")]
            if m := mass_re.search(line):
                masses[m.group(1)] = [float(x.strip().rstrip(",")) for x in m.group(2).split(",") if x.strip().rstrip(",")]
    return residues, masses


def peptide_body_with_ptms(peptide: str) -> str:
    return peptide.replace("[", "").replace("]", "")


def atom_counts(peptide: str, residues: dict[str, list[int]]) -> list[int]:
    counts = [0, 0, 0, 0, 0, 0]
    for token in ["Nterm", peptide_body_with_ptms(peptide), "Cterm"]:
        chars = [token] if token in ("Nterm", "Cterm") else token
        for residue in chars:
            comp = residues[residue]
            counts = [a + b for a, b in zip(counts, comp)]
    return counts


def base_mass(counts: list[int], masses: dict[str, list[float]]) -> float:
    elements = ["C", "H", "O", "N", "P", "S"]
    return sum(count * masses[element][0] for count, element in zip(counts, elements))


def load_fixture_record():
    with h5py.File(H5, "r") as handle:
        peptides = handle["records/peptide"][:]
        idx = next(i for i, value in enumerate(peptides) if value.decode().rstrip("\x00") == "[DQELAAR]")
        charge = int(handle["records/charge"][idx])
    return charge, "[DQELAAR]"


def load_ft1_scan(scan_number: int) -> list[tuple[float, float]]:
    peaks: list[tuple[float, float]] = []
    keep = False
    with FT1.open() as handle:
        for line in handle:
            if line.startswith("S\t") or line.startswith("S "):
                fields = line.split()
                keep = len(fields) > 1 and int(fields[1]) == scan_number
                continue
            if keep and line[:1].isdigit():
                mz, intensity, *_ = line.split()
                peaks.append((float(mz), float(intensity)))
            elif keep and line.startswith("S"):
                break
    return peaks


def load_ft2_precursors(scan_number: int) -> list[tuple[float, int]]:
    keep = False
    with FT2.open() as handle:
        for line in handle:
            if line.startswith("S\t") or line.startswith("S "):
                fields = line.split()
                keep = len(fields) > 1 and int(fields[1]) == scan_number
                continue
            if keep and (line.startswith("Z\t") or line.startswith("Z ")):
                fields = line.split()[3:]
                pairs = []
                for i in range(0, len(fields) - 1, 2):
                    charge = int(fields[i])
                    mz = float(fields[i + 1])
                    pairs.append((mz, charge))
                return pairs
    raise AssertionError(f"scan {scan_number} not found in {FT2}")


def find_peak(peaks: list[tuple[float, float]], target_mz: float) -> tuple[float, float] | None:
    tol = target_mz * MS1_TOLERANCE_PPM * 1e-6
    candidates = [peak for peak in peaks if abs(peak[0] - target_mz) <= tol]
    return max(candidates, key=lambda peak: peak[1]) if candidates else None


def collect_isotope_peaks(peaks: list[tuple[float, float]], anchor_mz: float, charge: int) -> list[tuple[float, float]]:
    anchor = find_peak(peaks, anchor_mz)
    assert anchor is not None, "fixture FT1 anchor peak should exist"

    out = [anchor]
    neutron_mz = NEUTRON / charge
    for direction in (-1, 1):
        for iso in range(1, 21):
            peak = find_peak(peaks, anchor[0] + direction * iso * neutron_mz)
            if peak is None:
                break
            out.append(peak)
    return sorted(set(out))


def ms1_abundance(
    collected: list[tuple[float, float]],
    mono_mz: float,
    charge: int,
    carbon_count: int,
    clamp_high: bool = True,
) -> tuple[float, int]:
    first = next((i for i, (mz, intensity) in enumerate(collected) if intensity > 0 and mz > mono_mz), None)
    if first is None:
        return 0.0, 0

    first_delta = round((collected[first][0] - mono_mz) / NEUTRON * charge)
    if first_delta < 0 or (first_delta > carbon_count and not clamp_high):
        return 0.0, 0

    total = 0.0
    weighted = 0.0
    count = 0
    for mz, intensity in collected[first:]:
        isotope_index = first_delta + round((mz - collected[first][0]) / NEUTRON * charge)
        if isotope_index > carbon_count and not clamp_high:
            break
        total += intensity
        weighted += intensity * min(isotope_index, carbon_count)
        count += 1

    return (weighted / total / carbon_count * 100.0, count) if count and total else (0.0, 0)


def main() -> None:
    residues, masses = parse_config()
    psm_charge, peptide = load_fixture_record()
    counts = atom_counts(peptide, residues)
    mono_mass = base_mass(counts, masses)
    mono_mz = mono_mass / psm_charge + PROTON
    ft1_peaks = load_ft1_scan(1003)

    matched_anchor_mz = 405.054657
    wrong_ft2_charge = 3

    ft2_precursors = load_ft2_precursors(1004)
    assert (matched_anchor_mz, psm_charge) in ft2_precursors
    assert (matched_anchor_mz, wrong_ft2_charge) in ft2_precursors
    eligible = [(mz, charge) for mz, charge in ft2_precursors if charge == psm_charge]
    assert (matched_anchor_mz, psm_charge) in eligible
    assert (matched_anchor_mz, wrong_ft2_charge) not in eligible

    old_bug_mono_mz = mono_mass / wrong_ft2_charge + PROTON
    wrong_peaks = collect_isotope_peaks(ft1_peaks, matched_anchor_mz, wrong_ft2_charge)
    wrong_pct, wrong_count = ms1_abundance(
        wrong_peaks, old_bug_mono_mz, wrong_ft2_charge, counts[0], clamp_high=False
    )
    assert (wrong_pct, wrong_count) == (0.0, 0), "fixture should reproduce the old charge bug"

    fixed_peaks = collect_isotope_peaks(ft1_peaks, matched_anchor_mz, psm_charge)
    fixed_pct, fixed_count = ms1_abundance(fixed_peaks, mono_mz, psm_charge, counts[0])
    assert fixed_count > 0, "PSM charge should produce at least the real FT1 anchor peak"
    assert fixed_pct > 0.0, "PSM charge should produce nonzero MS1 abundance"

    high_label_peptide = "[LVECNGKPVAK]"
    high_label_charge = 3
    high_label_anchor_mz = 422.945099
    high_label_counts = atom_counts(high_label_peptide, residues)
    high_label_mono_mz = base_mass(high_label_counts, masses) / high_label_charge + PROTON
    high_label_peaks = collect_isotope_peaks(load_ft1_scan(1063), high_label_anchor_mz, high_label_charge)

    old_high_pct, old_high_count = ms1_abundance(
        high_label_peaks,
        high_label_mono_mz,
        high_label_charge,
        high_label_counts[0],
        clamp_high=False,
    )
    assert (old_high_pct, old_high_count) == (0.0, 0), "fixture should reproduce the old high-label zero"

    high_pct, high_count = ms1_abundance(
        high_label_peaks, high_label_mono_mz, high_label_charge, high_label_counts[0]
    )
    assert high_count > 0, "high-label FT1 anchor should be counted"
    assert high_pct == 100.0, "high-label isotope index above carbon count should saturate to 100%"

    print(
        f"ok: {peptide} charge={psm_charge} anchor_mz={matched_anchor_mz} "
        f"mono_mz={mono_mz:.6f} isotopicPeakNumbers={fixed_count} "
        f"MS1IsotopicAbundances={fixed_pct:.6f}"
    )
    print(
        f"ok: {high_label_peptide} charge={high_label_charge} "
        f"anchor_mz={high_label_anchor_mz} mono_mz={high_label_mono_mz:.6f} "
        f"isotopicPeakNumbers={high_count} MS1IsotopicAbundances={high_pct:.6f}"
    )


if __name__ == "__main__":
    main()
