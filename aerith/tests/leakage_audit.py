#!/usr/bin/env python3
"""Black-box regression audit for Aerith held-out label isolation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import subprocess


FOLDS = 3


def spectrum_folds(rows: list[dict[str, str]]) -> list[int]:
    order = sorted(range(len(rows)), key=lambda i: (int(float(rows[i]["ScanNr"])), i))
    remaining = [0] * FOLDS
    left = len(order)
    for n in range(FOLDS, 0, -1):
        remaining[n - 1] = left // n
        left -= remaining[n - 1]
    seed = 1

    def draw() -> int:
        nonlocal seed
        seed = seed * 279470273 % 4294967291
        return seed % FOLDS

    folds = [0] * len(rows)
    begin = 0
    while begin < len(order):
        end = begin + 1
        scan = rows[order[begin]]["ScanNr"]
        while end < len(order) and rows[order[end]]["ScanNr"] == scan:
            end += 1
        fold = draw()
        while remaining[fold] <= 0:
            fold = draw()
        for position in range(begin, end):
            folds[order[position]] = fold
            remaining[fold] -= 1
        begin = end
    return folds


def run_aerith(binary: Path, pin: Path, prefix: Path) -> None:
    subprocess.run([
        str(binary), "--input", str(pin), "--output-prefix", str(prefix),
        "--q-threshold", "1",
    ], check=True, stdout=subprocess.DEVNULL)


def scores(prefix: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    for suffix in ("_target_psms.tsv", "_decoy_psms.tsv"):
        with Path(f"{prefix}{suffix}").open(newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                result[row["PSMId"]] = float(row["score"])
    return result


def rt_residuals(prefix: Path) -> dict[str, float]:
    with Path(f"{prefix}_filtered_psms.tsv").open(newline="") as stream:
        return {row["PSMId"]: float(row["sqrtAbsDeltaRT"])
                for row in csv.DictReader(stream, delimiter="\t")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aerith", required=True, type=Path)
    parser.add_argument("--pin", required=True, type=Path)
    parser.add_argument("--output", type=Path,
                        default=Path("/tmp/aerith-leakage-audit"))
    parser.add_argument("--fold", type=int, choices=range(FOLDS), default=0)
    parser.add_argument("--label-swaps", type=int, default=1000)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    args.aerith = args.aerith.resolve()
    args.pin = args.pin.resolve()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.pin.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    if fieldnames is None:
        raise SystemExit("PIN has no header")
    folds = spectrum_folds(rows)
    held_out_ids = {row["SpecId"] for row, fold in zip(rows, folds)
                    if fold == args.fold}
    held_out_targets = [i for i, fold in enumerate(folds)
                        if fold == args.fold and int(rows[i]["Label"]) == 1]
    held_out_decoys = [i for i, fold in enumerate(folds)
                       if fold == args.fold and int(rows[i]["Label"]) == -1]
    swaps = min(args.label_swaps, len(held_out_targets), len(held_out_decoys))
    changed_rows = set(held_out_targets[:swaps] + held_out_decoys[:swaps])

    perturbed = args.output / "heldout_labels_flipped.pin"
    with perturbed.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows):
            changed = dict(row)
            if index in changed_rows:
                changed["Label"] = "-1" if int(row["Label"]) == 1 else "1"
            writer.writerow(changed)

    original_prefix = args.output / "original"
    changed_prefix = args.output / "perturbed"
    run_aerith(args.aerith, args.pin, original_prefix)
    run_aerith(args.aerith, perturbed, changed_prefix)

    original_scores, changed_scores = scores(original_prefix), scores(changed_prefix)
    original_residuals = rt_residuals(original_prefix)
    changed_residuals = rt_residuals(changed_prefix)
    score_delta = max(abs(original_scores[i] - changed_scores[i]) for i in held_out_ids)
    rt_delta = max(abs(original_residuals[i] - changed_residuals[i]) for i in held_out_ids)
    passed = score_delta <= args.tolerance and rt_delta <= args.tolerance
    report = "\n".join([
        "Aerith held-out label leakage audit",
        "===================================",
        f"Fold:                              {args.fold}",
        f"Held-out PSMs:                     {len(held_out_ids)}",
        f"Changed held-out labels:           {len(changed_rows)}",
        f"Maximum held-out score delta:      {score_delta:.12g}",
        f"Maximum held-out RT feature delta: {rt_delta:.12g}",
        f"Tolerance:                         {args.tolerance:.12g}",
        f"Result:                            {'PASS' if passed else 'FAIL'}",
        "",
    ])
    (args.output / "audit.txt").write_text(report)
    print(report, end="")
    if not passed:
        raise SystemExit("held-out label leakage detected")


if __name__ == "__main__":
    main()
