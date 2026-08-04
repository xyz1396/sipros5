import argparse
import csv
import ctypes
import math
import mmap
import random
import re
import struct
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import h5py
import numpy as np


MAX_PLOT_RECORDS = 8


@dataclass
class SpectraRecord:
    source_file: str
    record_kind: str
    target_sip_abundance_pct: float | None
    sip_label: str
    psm_id: str
    retention: str
    charge: int
    peptide: str
    proteins: str
    precursor_mz: np.ndarray
    precursor_intensity: np.ndarray
    fragment_mz: np.ndarray
    theoretical_intensity: np.ndarray
    experimental_intensity: np.ndarray
    ion_kinds: list[str]
    ion_positions: np.ndarray
    source_index: int
    predicted_fragment_mz: np.ndarray | None = None
    predicted_fragment_intensity: np.ndarray | None = None
    predicted_ion_kinds: list[str] | None = None
    predicted_ion_positions: np.ndarray | None = None
    predicted_ion_charges: np.ndarray | None = None
    predicted_rt: float | None = None
    raw_fragment_mz: np.ndarray | None = None
    raw_fragment_intensity: np.ndarray | None = None


@dataclass
class PredictedSpectrum:
    mz: np.ndarray
    intensity: np.ndarray
    ion_kinds: list[str]
    ion_positions: np.ndarray
    ion_charges: np.ndarray


SPECTRUM_CACHE_MAGIC = 0x4145525350454332
RT_CACHE_MAGIC = 0x4145525254505232


class SfiHeader(ctypes.Structure):
    _fields_ = [
        ("magic", ctypes.c_char * 8),
        ("version", ctypes.c_uint32),
        ("endian", ctypes.c_uint32),
        ("headerSize", ctypes.c_uint32),
        ("recordSize", ctypes.c_uint32),
        ("precursorPeakSize", ctypes.c_uint32),
        ("fragmentPeakSize", ctypes.c_uint32),
        ("productPostingSize", ctypes.c_uint32),
        ("rtBinSize", ctypes.c_uint32),
        ("recordsPerBlock", ctypes.c_uint32),
        ("blockCount", ctypes.c_uint32),
        ("sipIsotopeMassNumber", ctypes.c_uint32),
        ("label", ctypes.c_int32),
        ("sipAtom", ctypes.c_char),
        ("reservedChars", ctypes.c_char * 7),
        ("fileSize", ctypes.c_uint64),
        ("recordCount", ctypes.c_uint64),
        ("precursorCount", ctypes.c_uint64),
        ("fragmentCount", ctypes.c_uint64),
        ("rtBinCount", ctypes.c_uint64),
        ("productPostingCount", ctypes.c_uint64),
        ("stringBytes", ctypes.c_uint64),
        ("recordOffset", ctypes.c_uint64),
        ("precursorOffset", ctypes.c_uint64),
        ("fragmentOffset", ctypes.c_uint64),
        ("blockRtBinOffset", ctypes.c_uint64),
        ("blockProductOffset", ctypes.c_uint64),
        ("rtBinOffset", ctypes.c_uint64),
        ("productPostingOffset", ctypes.c_uint64),
        ("stringOffset", ctypes.c_uint64),
        ("fragmentBinWidth", ctypes.c_double),
        ("rtBinWidthMinutes", ctypes.c_double),
        ("targetSipAbundancePct", ctypes.c_double),
        ("probabilityCutoff", ctypes.c_double),
        ("generationPpmTolerance", ctypes.c_double),
        ("minimumMatchedEnvelopes", ctypes.c_uint64),
        ("payloadChecksum", ctypes.c_uint64),
        ("chemistryProfileId", ctypes.c_char * 96),
        ("recordKind", ctypes.c_char * 24),
        ("envelopeTopN", ctypes.c_uint32),
        ("reservedEnvelope", ctypes.c_uint32),
        ("reserved", ctypes.c_uint64 * 7),
    ]


class SfiRecord(ctypes.Structure):
    _fields_ = [
        ("topPrecursorMz", ctypes.c_double),
        ("topPrecursorIntensity", ctypes.c_double),
        ("sumPrecursorIntensity", ctypes.c_double),
        ("retentionMinutes", ctypes.c_double),
        ("sipAbundancePct", ctypes.c_double),
        ("precursorOffset", ctypes.c_uint64),
        ("fragmentOffset", ctypes.c_uint64),
        ("psmIdOffset", ctypes.c_uint64),
        ("peptideOffset", ctypes.c_uint64),
        ("proteinsOffset", ctypes.c_uint64),
        ("precursorCount", ctypes.c_uint32),
        ("fragmentCount", ctypes.c_uint32),
        ("psmIdSize", ctypes.c_uint32),
        ("peptideSize", ctypes.c_uint32),
        ("proteinsSize", ctypes.c_uint32),
        ("generationOrdinal", ctypes.c_uint32),
        ("charge", ctypes.c_int32),
        ("reserved", ctypes.c_uint32),
    ]


PRECURSOR_DTYPE = np.dtype([("mz", "<f8"), ("intensity", "<f8")])
FRAGMENT_DTYPE = np.dtype({
    "names": ["mz_bin", "theory", "experiment", "position", "kind", "reserved"],
    "formats": ["<u4", "<f4", "<f4", "<u2", "u1", "u1"],
    "offsets": [0, 4, 8, 12, 14, 15],
    "itemsize": 16,
})


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Visualize a Sipros SFI library together with its mandatory Aerith "
            "predicted-spectrum and RT caches and mandatory real-sample Raxport HDF5. "
            "The script never substitutes a theoretical-spectra HDF5 file for SFI or "
            "reruns prediction when a cache entry is missing."
        )
    )
    parser.add_argument(
        "--sfi", "--input", dest="sfi", required=True,
        help="SFI file, or a directory containing target/decoy .sfi files.",
    )
    parser.add_argument(
        "--spectrum-cache", required=True,
        help="Aerith regular_search_predictions.spectrum file.",
    )
    parser.add_argument(
        "--rt-cache", required=True,
        help="Aerith regular_search_predictions.rt file.",
    )
    parser.add_argument(
        "--sample-hdf5", required=True,
        help=("Real-sample Raxport HDF5 file, or a directory recursively containing "
              "one <sample>.h5 file per SFI source sample."),
    )
    parser.add_argument(
        "--output",
        default="test/matched_experimental_spectra_top.pdf",
        help="Output plot path.",
    )
    parser.add_argument(
        "--tsv-output",
        default="",
        help="Output TSV path for selected peak rows. Defaults to the plot path with .tsv suffix.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of records to select; plots at most 8.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=1000,
        help="Number of records to sample before ranking. Use 0 to rank all records.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed used when sampling records.",
    )
    parser.add_argument(
        "--rank-by",
        choices=["matched-envelopes", "matched-peaks", "source-order", "random"],
        default="matched-envelopes",
        help="How selected records are ranked after sampling.",
    )
    parser.add_argument(
        "--include-low-matches",
        action="store_true",
        help=(
            "Split --top-k between the highest- and lowest-match records "
            "instead of plotting only the highest matches."
        ),
    )
    parser.add_argument(
        "--psm-id",
        action="append",
        default=[],
        help="Plot a specific PSM id. Can be supplied multiple times; bypasses sampling/ranking.",
    )
    parser.add_argument(
        "--min-matched-envelopes",
        type=int,
        default=0,
        help="Only consider records with at least this many matched fragment envelopes.",
    )
    parser.add_argument(
        "--mz-min",
        type=float,
        default=None,
        help="Optional lower m/z bound for the fragment panel.",
    )
    parser.add_argument(
        "--mz-max",
        type=float,
        default=None,
        help="Optional upper m/z bound for the fragment panel.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure after saving it.",
    )
    return parser.parse_args()


def _decode_fixed(value):
    return bytes(value).split(b"\0", 1)[0].decode("utf-8")


class SfiIndex:
    def __init__(self, path):
        self.path = Path(path)
        self._stream = self.path.open("rb")
        self._mapping = mmap.mmap(self._stream.fileno(), 0, access=mmap.ACCESS_READ)
        if len(self._mapping) < ctypes.sizeof(SfiHeader):
            raise ValueError(f"{self.path}: truncated SFI header")
        self.header = SfiHeader.from_buffer_copy(
            self._mapping[: ctypes.sizeof(SfiHeader)]
        )
        if bytes(self.header.magic) != b"SIPSFI05":
            raise ValueError(
                f"{self.path}: not an SFI v5 file; theoretical-spectra HDF5 fallback is disabled"
            )
        if self.header.version != 5 or self.header.endian != 0x01020304:
            raise ValueError(f"{self.path}: unsupported SFI version or byte order")
        expected = {
            "header": ctypes.sizeof(SfiHeader),
            "record": ctypes.sizeof(SfiRecord),
            "precursor": PRECURSOR_DTYPE.itemsize,
            "fragment": FRAGMENT_DTYPE.itemsize,
        }
        observed = {
            "header": self.header.headerSize,
            "record": self.header.recordSize,
            "precursor": self.header.precursorPeakSize,
            "fragment": self.header.fragmentPeakSize,
        }
        if expected != observed or self.header.fileSize != len(self._mapping):
            raise ValueError(f"{self.path}: unsupported or inconsistent SFI layout")
        self.record_kind = _decode_fixed(self.header.recordKind)
        atom = self.header.sipAtom.decode("ascii") if self.header.sipAtom else ""
        self.sip_label = (
            f"{atom}{self.header.sipIsotopeMassNumber}" if atom else ""
        )

    def close(self):
        self._mapping.close()
        self._stream.close()

    @property
    def record_count(self):
        return int(self.header.recordCount)

    def _record(self, index):
        if index < 0 or index >= self.record_count:
            raise IndexError(index)
        offset = self.header.recordOffset + index * ctypes.sizeof(SfiRecord)
        return SfiRecord.from_buffer_copy(
            self._mapping[offset : offset + ctypes.sizeof(SfiRecord)]
        )

    def _text(self, offset, size):
        begin = self.header.stringOffset + int(offset)
        end = begin + int(size)
        if end > self.header.stringOffset + self.header.stringBytes:
            raise ValueError(f"{self.path}: SFI string lies outside its table")
        return self._mapping[begin:end].decode("utf-8")

    def psm_id(self, index):
        record = self._record(index)
        return self._text(record.psmIdOffset, record.psmIdSize)

    def materialize(self, index):
        record = self._record(index)
        precursor_offset = (
            self.header.precursorOffset +
            int(record.precursorOffset) * PRECURSOR_DTYPE.itemsize
        )
        fragment_offset = (
            self.header.fragmentOffset +
            int(record.fragmentOffset) * FRAGMENT_DTYPE.itemsize
        )
        precursors = np.frombuffer(
            self._mapping, dtype=PRECURSOR_DTYPE,
            count=int(record.precursorCount), offset=precursor_offset,
        ).copy()
        fragments = np.frombuffer(
            self._mapping, dtype=FRAGMENT_DTYPE,
            count=int(record.fragmentCount), offset=fragment_offset,
        ).copy()
        return SpectraRecord(
            source_file=self.path.name,
            record_kind=self.record_kind,
            target_sip_abundance_pct=float(record.sipAbundancePct),
            sip_label=self.sip_label,
            psm_id=self._text(record.psmIdOffset, record.psmIdSize),
            retention=f"{float(record.retentionMinutes):.10g}",
            charge=int(record.charge),
            peptide=self._text(record.peptideOffset, record.peptideSize),
            proteins=self._text(record.proteinsOffset, record.proteinsSize),
            precursor_mz=precursors["mz"].astype(float),
            precursor_intensity=precursors["intensity"].astype(float),
            fragment_mz=fragments["mz_bin"].astype(float) * 0.001,
            theoretical_intensity=fragments["theory"].astype(float),
            experimental_intensity=fragments["experiment"].astype(float),
            ion_kinds=[chr(int(value)) for value in fragments["kind"]],
            ion_positions=fragments["position"].astype(int),
            source_index=index,
        )


def _sfi_paths(path):
    path = Path(path)
    if path.is_dir():
        files = sorted(item for item in path.iterdir()
                       if item.is_file() and item.suffix.lower() == ".sfi")
        if not files:
            raise ValueError(f"{path}: no .sfi files found")
        return files
    if path.suffix.lower() != ".sfi":
        raise ValueError(
            f"{path}: expected SFI input; theoretical-spectra HDF5 fallback is disabled"
        )
    return [path]


def load_sfi_candidates(path, args):
    indexes = [SfiIndex(item) for item in _sfi_paths(path)]
    try:
        if args.psm_id:
            requested = set(args.psm_id)
            records = []
            for index in indexes:
                for row in range(index.record_count):
                    if index.psm_id(row) in requested:
                        records.append(index.materialize(row))
            return records, sum(index.record_count for index in indexes), len(indexes)

        references = [
            (library, row)
            for library, index in enumerate(indexes)
            for row in range(index.record_count)
        ] if args.sample_size == 0 else None
        total = sum(index.record_count for index in indexes)
        if references is None:
            count = min(args.sample_size, total)
            selected = sorted(random.Random(args.seed).sample(range(total), count))
            references = []
            library = 0
            begin = 0
            for global_row in selected:
                while global_row >= begin + indexes[library].record_count:
                    begin += indexes[library].record_count
                    library += 1
                references.append((library, global_row - begin))
        return ([indexes[library].materialize(row)
                 for library, row in references], total, len(indexes))
    finally:
        for index in indexes:
            index.close()


def _peptide_body(peptide):
    first = peptide.find("[")
    last = peptide.rfind("]")
    if (first >= 0 and last > first and peptide.find("[", first + 1) < 0
            and first <= 2 and len(peptide) - last <= 3):
        return peptide[first + 1:last]
    first = peptide.find(".")
    last = peptide.rfind(".")
    if first >= 0 and last > first:
        return peptide[first + 1:last]
    return peptide


def _spectrum_key(record):
    return (_peptide_body(record.peptide) + "\x1f" + str(record.charge)).encode()


def _numeric_modification(body, position):
    if position >= len(body) or body[position] != "[":
        return None, position
    close = body.find("]", position + 1)
    if close < 0:
        return None, position
    try:
        return float(body[position + 1:close]), close + 1
    except ValueError:
        return None, position


def _rt_token_key(peptide):
    body = _peptide_body(peptide)
    residue_tokens = {
        "G": 3, "A": 4, "V": 5, "I": 6, "L": 7, "P": 8,
        "F": 9, "W": 10, "M": 11, "X": 11, "S": 13, "T": 14,
        "Y": 15, "Q": 16, "E": 17, "N": 18, "D": 19, "K": 20,
        "O": 20, "R": 21, "H": 22, "C": 24, "U": 24,
    }
    modification_symbols = set("~!@><%^&*()/$")
    position = 0
    acetylated = body.startswith("%")
    if acetylated:
        position = 1
    else:
        shift, next_position = _numeric_modification(body, position)
        if shift is not None and abs(shift - 42.0106) < 0.02:
            acetylated = True
            position = next_position
    tokens = [29 if acetylated else 1]
    residues = 0
    while position < len(body):
        residue = body[position]
        position += 1
        if not ("A" <= residue <= "Z"):
            continue
        if residue not in residue_tokens:
            raise ValueError(f"Unsupported DIA-NN RT residue in {peptide}: {residue}")
        modification = ""
        if position < len(body) and body[position] in modification_symbols:
            modification = body[position]
            position += 1
        numeric, next_position = _numeric_modification(body, position)
        if numeric is not None:
            position = next_position
        token = residue_tokens[residue]
        if residue == "C" and modification != "(" and (
                numeric is None or abs(numeric - 57.0215) < 0.02):
            token = 25
        elif residue == "M" and (
                modification == "~" or
                (numeric is not None and abs(numeric - 15.9949) < 0.02)):
            token = 26
        elif residue in "STY" and (
                modification in "@><" or
                (numeric is not None and abs(numeric - 79.9663) < 0.02)):
            token = {"S": 31, "T": 32, "Y": 33}[residue]
        elif residue == "K" and (
                modification == "&" or
                (numeric is not None and abs(numeric - 28.0313) < 0.02)):
            token = 39
        tokens.append(token)
        residues += 1
    if residues < 5:
        raise ValueError(f"DIA-NN RT prediction requires five residues: {peptide}")
    tokens.append(2)
    return bytes(tokens)


def _read_cache_header(stream, expected_magic, path):
    header = stream.read(24)
    if len(header) != 24:
        raise ValueError(f"{path}: truncated prediction-cache header")
    magic, _fingerprint, count = struct.unpack("<QQQ", header)
    if magic != expected_magic:
        raise ValueError(f"{path}: incompatible prediction-cache format")
    if count > 1_000_000_000:
        raise ValueError(f"{path}: unreasonable prediction-cache entry count")
    return count


def load_spectrum_cache(path, requested):
    result = {}
    with Path(path).open("rb") as stream:
        count = _read_cache_header(stream, SPECTRUM_CACHE_MAGIC, path)
        fragment_struct = struct.Struct("<ffcIi")
        for _ in range(count):
            raw = stream.read(4)
            if len(raw) != 4:
                raise ValueError(f"{path}: truncated spectrum-cache key size")
            key_size, = struct.unpack("<I", raw)
            if key_size > 65536:
                raise ValueError(f"{path}: invalid spectrum-cache key size")
            key = stream.read(key_size)
            raw = stream.read(4)
            if len(key) != key_size or len(raw) != 4:
                raise ValueError(f"{path}: truncated spectrum-cache entry")
            fragment_count, = struct.unpack("<I", raw)
            if fragment_count > 10000:
                raise ValueError(f"{path}: invalid spectrum-cache fragment count")
            if key not in requested:
                stream.seek(fragment_count * fragment_struct.size, 1)
                continue
            fragments = []
            for _fragment in range(fragment_count):
                raw = stream.read(fragment_struct.size)
                if len(raw) != fragment_struct.size:
                    raise ValueError(f"{path}: truncated predicted fragment")
                fragments.append(fragment_struct.unpack(raw))
            result[key] = PredictedSpectrum(
                mz=np.asarray([value[0] for value in fragments], dtype=float),
                intensity=np.asarray([value[1] for value in fragments], dtype=float),
                ion_kinds=[value[2].decode("ascii") for value in fragments],
                ion_positions=np.asarray([value[3] for value in fragments], dtype=int),
                ion_charges=np.asarray([value[4] for value in fragments], dtype=int),
            )
    missing = requested.difference(result)
    if missing:
        raise ValueError(
            f"{path}: {len(missing)} selected SFI peptide-charge forms are absent; "
            "prediction fallback is disabled"
        )
    return result


def load_rt_cache(path, requested):
    result = {}
    with Path(path).open("rb") as stream:
        count = _read_cache_header(stream, RT_CACHE_MAGIC, path)
        for _ in range(count):
            raw = stream.read(4)
            if len(raw) != 4:
                raise ValueError(f"{path}: truncated RT-cache key size")
            key_size, = struct.unpack("<I", raw)
            if key_size > 65536:
                raise ValueError(f"{path}: invalid RT-cache key size")
            key = stream.read(key_size)
            raw = stream.read(4)
            if len(key) != key_size or len(raw) != 4:
                raise ValueError(f"{path}: truncated RT-cache entry")
            if key in requested:
                result[key] = struct.unpack("<f", raw)[0]
    missing = requested.difference(result)
    if missing:
        raise ValueError(
            f"{path}: {len(missing)} selected SFI peptides lack predicted RT; "
            "prediction fallback is disabled"
        )
    return result


def attach_predictions(records, spectrum_cache_path, rt_cache_path):
    spectrum_keys = {_spectrum_key(record) for record in records}
    rt_keys = {_rt_token_key(record.peptide) for record in records}
    spectra = load_spectrum_cache(spectrum_cache_path, spectrum_keys)
    retention = load_rt_cache(rt_cache_path, rt_keys)
    for record in records:
        predicted = spectra[_spectrum_key(record)]
        record.predicted_fragment_mz = predicted.mz
        record.predicted_fragment_intensity = predicted.intensity
        record.predicted_ion_kinds = predicted.ion_kinds
        record.predicted_ion_positions = predicted.ion_positions
        record.predicted_ion_charges = predicted.ion_charges
        record.predicted_rt = retention[_rt_token_key(record.peptide)]


def _sample_name_and_scan(psm_id):
    source_id = re.sub(r"_\d{3}\.\d{3}Pct$", "", psm_id.strip())
    if source_id.startswith("DECOY_"):
        source_id = source_id[len("DECOY_"):]
    parts = source_id.split(".")
    if len(parts) < 3:
        raise ValueError(f"Cannot parse sample/scan from SFI PSM id: {psm_id}")
    for value in reversed(parts[:-1]):
        try:
            scan = int(value)
        except ValueError:
            continue
        if scan > 0:
            return parts[0], scan
    raise ValueError(f"Cannot parse scan from SFI PSM id: {psm_id}")


def _sample_hdf5_paths(path):
    path = Path(path)
    files = (sorted(
        list(path.rglob("*.h5")) + list(path.rglob("*.hdf5"))
    ) if path.is_dir() else [path])
    if not files or any(item.suffix.lower() not in {".h5", ".hdf5"} for item in files):
        raise ValueError(f"{path}: expected real-sample Raxport HDF5 input")
    by_sample = {}
    for item in files:
        if item.stem in by_sample:
            raise ValueError(f"Duplicate real-sample HDF5 basename: {item.stem}")
        by_sample[item.stem] = item
    return by_sample


def attach_raw_spectra(records, sample_hdf5_path):
    paths = _sample_hdf5_paths(sample_hdf5_path)
    requested = {}
    for record in records:
        sample, scan = _sample_name_and_scan(record.psm_id)
        if sample not in paths:
            raise ValueError(f"No real-sample HDF5 for SFI sample {sample}")
        requested.setdefault(sample, {}).setdefault(scan, []).append(record)
    for sample, scans in requested.items():
        with h5py.File(paths[sample], "r") as handle:
            scan_numbers = np.asarray(handle["scans/scan_number"][:], dtype=int)
            ms_order = np.asarray(handle["scans/ms_order"][:], dtype=int)
            peak_start = np.asarray(handle["scans/peak_start"][:], dtype=np.int64)
            peak_count = np.asarray(handle["scans/peak_count"][:], dtype=int)
            row_by_scan = {
                int(scan): row for row, scan in enumerate(scan_numbers)
                if int(ms_order[row]) == 2
            }
            for scan, scan_records in scans.items():
                if scan not in row_by_scan:
                    raise ValueError(f"{paths[sample]}: MS2 scan {scan} is absent")
                row = row_by_scan[scan]
                start = int(peak_start[row])
                count = int(peak_count[row])
                mz = np.asarray(handle["peaks/mz"][start:start + count], dtype=float)
                intensity = np.asarray(
                    handle["peaks/intensity"][start:start + count], dtype=float
                )
                maximum = float(np.max(intensity)) if intensity.size else 0.0
                if maximum > 0.0:
                    intensity = intensity / maximum
                for record in scan_records:
                    record.raw_fragment_mz = mz.copy()
                    record.raw_fragment_intensity = intensity.copy()


def matched_peak_count(record):
    return int(np.count_nonzero(record.experimental_intensity > 0.0))


def envelope_keys(record):
    return list(zip(record.ion_kinds, record.ion_positions.tolist()))


def matched_envelope_count(record):
    matched = set()
    for key, intensity in zip(envelope_keys(record), record.experimental_intensity):
        if intensity > 0.0:
            matched.add(key)
    return len(matched)


def retained_envelope_count(record):
    return len(set(envelope_keys(record)))


def _select_records(records, args):
    if args.psm_id:
        by_id = {record.psm_id: record for record in records}
        missing = [psm_id for psm_id in args.psm_id if psm_id not in by_id]
        if missing:
            raise SystemExit(f"PSM id not found: {', '.join(missing)}")
        return [by_id[psm_id] for psm_id in args.psm_id]

    filtered = [
        record
        for record in records
        if matched_envelope_count(record) >= args.min_matched_envelopes
    ]
    sampled = list(filtered)
    if args.rank_by == "random":
        rng = random.Random(args.seed)
        rng.shuffle(sampled)
        return sampled[: args.top_k]
    if args.rank_by == "source-order":
        ranked = sorted(sampled, key=lambda record: record.source_index)
    elif args.rank_by == "matched-peaks":
        ranked = sorted(
            sampled,
            key=lambda record: (matched_peak_count(record), matched_envelope_count(record), -record.source_index),
            reverse=True,
        )
    else:
        ranked = sorted(
            sampled,
            key=lambda record: (matched_envelope_count(record), matched_peak_count(record), -record.source_index),
            reverse=True,
        )
    if args.include_low_matches and args.rank_by in {"matched-envelopes", "matched-peaks"}:
        high_count = (args.top_k + 1) // 2
        low_count = args.top_k - high_count
        high = ranked[:high_count]
        if args.rank_by == "matched-peaks":
            low_ranked = sorted(
                sampled,
                key=lambda record: (
                    matched_peak_count(record),
                    matched_envelope_count(record),
                    record.source_index,
                ),
            )
        else:
            low_ranked = sorted(
                sampled,
                key=lambda record: (
                    matched_envelope_count(record),
                    matched_peak_count(record),
                    record.source_index,
                ),
            )
        selected_indices = {
            (record.source_file, record.source_index) for record in high
        }
        low = [
            record
            for record in low_ranked
            if (record.source_file, record.source_index) not in selected_indices
        ][:low_count]
        return high + low
    return ranked[: args.top_k]


def _format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.10g}"
        return ""
    return str(value)


def _safe_float(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _truncate_text(text, max_len=48):
    text = str(text or "")
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def _plot_peak_series(ax, mz, intensity, color, label, direction=1.0, linewidth=0.5, alpha=0.72):
    if mz.size == 0:
        return
    ax.vlines(
        mz,
        0.0,
        direction * intensity,
        color=color,
        alpha=alpha,
        linewidth=linewidth,
        label=label,
    )


def _set_mz_xlim(ax, mz_values, mz_min=None, mz_max=None):
    if mz_min is not None or mz_max is not None:
        if mz_min is None:
            mz_min = float(np.min(mz_values)) if mz_values.size else 0.0
        if mz_max is None:
            mz_max = float(np.max(mz_values)) if mz_values.size else mz_min + 1.0
        ax.set_xlim(mz_min, mz_max)
        return
    if mz_values.size == 0:
        ax.set_xlim(0.0, 1.0)
        return
    min_mz = float(np.min(mz_values))
    max_mz = float(np.max(mz_values))
    if math.isclose(min_mz, max_mz):
        pad = max(1.0, 0.02 * max_mz)
    else:
        pad = max(1.0, 0.02 * (max_mz - min_mz))
    ax.set_xlim(min_mz - pad, max_mz + pad)


def _dedup_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    dedup = {}
    for handle, label in zip(handles, labels):
        dedup.setdefault(label, handle)
    ax.legend(dedup.values(), dedup.keys(), loc="upper right", fontsize=8)


def _plot_precursor(ax, record):
    max_intensity = float(record.precursor_intensity.max()) if record.precursor_intensity.size else 1.0
    _plot_peak_series(
        ax,
        record.precursor_mz,
        record.precursor_intensity,
        color="tab:green",
        label="Theoretical precursor",
        linewidth=1.2,
        alpha=0.85,
    )
    ax.set_ylim(0.0, max(1.05, max_intensity * 1.08))
    _set_mz_xlim(ax, record.precursor_mz)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity")
    ax.set_title("Precursor Envelope")
    _dedup_legend(ax)


def _plot_fragments(ax, record, mz_min=None, mz_max=None):
    if mz_min is None and mz_max is None:
        mask = np.ones(record.fragment_mz.shape[0], dtype=bool)
    else:
        mask = np.ones(record.fragment_mz.shape[0], dtype=bool)
        if mz_min is not None:
            mask &= record.fragment_mz >= mz_min
        if mz_max is not None:
            mask &= record.fragment_mz <= mz_max

    mz = record.fragment_mz[mask]
    theory = record.theoretical_intensity[mask]
    experiment = record.experimental_intensity[mask]
    kinds = np.asarray(record.ion_kinds, dtype=object)[mask]

    matched_mask = experiment > 0.0
    b_mask = kinds == "b"
    y_mask = kinds == "y"

    raw_mz = np.asarray(record.raw_fragment_mz, dtype=float)
    raw_intensity = np.asarray(record.raw_fragment_intensity, dtype=float)
    raw_mask = np.ones(raw_mz.shape[0], dtype=bool)
    if mz_min is not None:
        raw_mask &= raw_mz >= mz_min
    if mz_max is not None:
        raw_mask &= raw_mz <= mz_max
    _plot_peak_series(
        ax, raw_mz[raw_mask], raw_intensity[raw_mask], color="0.72",
        label="Real sample MS/MS", linewidth=0.35, alpha=0.5,
    )

    if record.record_kind.lower() == "theoretical":
        _plot_peak_series(
            ax,
            mz[b_mask],
            theory[b_mask],
            color="tab:blue",
            label="Theoretical b",
            linewidth=0.45,
            alpha=0.7,
        )
        _plot_peak_series(
            ax,
            mz[y_mask],
            theory[y_mask],
            color="tab:purple",
            label="Theoretical y",
            linewidth=0.45,
            alpha=0.7,
        )
        max_intensity = float(theory.max()) if theory.size else 1.0
        ax.set_ylim(0.0, max(1.05, max_intensity * 1.08))
        _set_mz_xlim(ax, mz, mz_min, mz_max)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Relative intensity")
        ax.set_title("Theoretical Fragments")
        _dedup_legend(ax)
        return

    _plot_peak_series(
        ax,
        mz[matched_mask & b_mask],
        experiment[matched_mask & b_mask],
        color="tab:blue",
        label="Observed-apex-supported b model",
        linewidth=0.5,
    )
    _plot_peak_series(
        ax,
        mz[matched_mask & y_mask],
        experiment[matched_mask & y_mask],
        color="tab:purple",
        label="Observed-apex-supported y model",
        linewidth=0.5,
    )
    _plot_peak_series(
        ax,
        mz[~matched_mask],
        theory[~matched_mask],
        color="0.78",
        label="Unmatched theoretical",
        direction=-1.0,
        linewidth=0.35,
        alpha=0.55,
    )
    _plot_peak_series(
        ax,
        mz[matched_mask],
        theory[matched_mask],
        color="tab:orange",
        label="Matched theoretical",
        direction=-1.0,
        linewidth=0.4,
        alpha=0.65,
    )

    ax.axhline(0.0, color="0.5", linewidth=0.8)
    ax.set_ylim(-1.05, 1.05)
    _set_mz_xlim(ax, mz, mz_min, mz_max)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity")
    ax.set_title("SFI Envelope Model vs Real Sample MS/MS")
    ax.text(0.01, 0.96, "Real MS/MS + SFI apex model", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    ax.text(0.01, 0.04, "Theoretical", transform=ax.transAxes, ha="left", va="bottom", fontsize=8)
    _dedup_legend(ax)


def _plot_cached_prediction(ax, record, mz_min=None, mz_max=None):
    mz = np.asarray(record.predicted_fragment_mz, dtype=float)
    intensity = np.asarray(record.predicted_fragment_intensity, dtype=float)
    kinds = np.asarray(record.predicted_ion_kinds, dtype=object)
    mask = np.ones(mz.shape[0], dtype=bool)
    if mz_min is not None:
        mask &= mz >= mz_min
    if mz_max is not None:
        mask &= mz <= mz_max
    b_mask = kinds == "b"
    y_mask = kinds == "y"
    _plot_peak_series(
        ax, mz[mask & b_mask], intensity[mask & b_mask], "tab:blue",
        "Cached DIA-NN b", linewidth=0.7, alpha=0.8,
    )
    _plot_peak_series(
        ax, mz[mask & y_mask], intensity[mask & y_mask], "tab:purple",
        "Cached DIA-NN y", linewidth=0.7, alpha=0.8,
    )
    raw_mz = np.asarray(record.raw_fragment_mz, dtype=float)
    raw_intensity = np.asarray(record.raw_fragment_intensity, dtype=float)
    raw_mask = np.ones(raw_mz.shape[0], dtype=bool)
    if mz_min is not None:
        raw_mask &= raw_mz >= mz_min
    if mz_max is not None:
        raw_mask &= raw_mz <= mz_max
    _plot_peak_series(
        ax, raw_mz[raw_mask], raw_intensity[raw_mask], "0.55",
        "Real sample MS/MS", direction=-1.0, linewidth=0.35, alpha=0.45,
    )
    ax.axhline(0.0, color="0.5", linewidth=0.8)
    ax.set_ylim(-1.05, 1.05)
    _set_mz_xlim(ax, np.concatenate((mz[mask], raw_mz[raw_mask])), mz_min, mz_max)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity")
    ax.set_title("Cached DIA-NN Prediction vs Real Sample MS/MS")
    _dedup_legend(ax)


def _write_selected_tsv(records, output_path):
    fieldnames = [
        "record_rank",
        "source_file",
        "record_kind",
        "sip_label",
        "target_sip_abundance_pct",
        "psm_id",
        "retention",
        "cached_predicted_rt",
        "charge",
        "peptide",
        "proteins",
        "matched_envelopes",
        "retained_envelopes",
        "supported_modeled_peaks",
        "peak_index",
        "peak_kind",
        "ion_kind",
        "ion_position",
        "fragment_charge",
        "mz",
        "theoretical_intensity",
        "experimental_intensity",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record_rank, record in enumerate(records, start=1):
            base = {
                "record_rank": record_rank,
                "source_file": record.source_file,
                "record_kind": record.record_kind,
                "sip_label": record.sip_label,
                "target_sip_abundance_pct": _format_value(record.target_sip_abundance_pct),
                "psm_id": record.psm_id,
                "retention": record.retention,
                "cached_predicted_rt": _format_value(record.predicted_rt),
                "charge": record.charge,
                "peptide": record.peptide,
                "proteins": record.proteins,
                "matched_envelopes": matched_envelope_count(record),
                "retained_envelopes": retained_envelope_count(record),
                "supported_modeled_peaks": matched_peak_count(record),
            }
            for peak_index, mz in enumerate(record.precursor_mz):
                row = dict(base)
                row.update(
                    {
                        "peak_index": peak_index,
                        "peak_kind": "precursor",
                        "mz": _format_value(float(mz)),
                        "theoretical_intensity": _format_value(float(record.precursor_intensity[peak_index])),
                    }
                )
                writer.writerow(row)
                row_count += 1
            for peak_index, mz in enumerate(record.predicted_fragment_mz):
                row = dict(base)
                row.update(
                    {
                        "peak_index": peak_index,
                        "peak_kind": "cached_predicted_fragment",
                        "ion_kind": record.predicted_ion_kinds[peak_index],
                        "ion_position": int(record.predicted_ion_positions[peak_index]),
                        "fragment_charge": int(record.predicted_ion_charges[peak_index]),
                        "mz": _format_value(float(mz)),
                        "theoretical_intensity": _format_value(
                            float(record.predicted_fragment_intensity[peak_index])
                        ),
                    }
                )
                writer.writerow(row)
                row_count += 1
            for peak_index, mz in enumerate(record.raw_fragment_mz):
                row = dict(base)
                row.update(
                    {
                        "peak_index": peak_index,
                        "peak_kind": "real_sample_fragment",
                        "mz": _format_value(float(mz)),
                        "experimental_intensity": _format_value(
                            float(record.raw_fragment_intensity[peak_index])
                        ),
                    }
                )
                writer.writerow(row)
                row_count += 1
            for peak_index, mz in enumerate(record.fragment_mz):
                row = dict(base)
                row.update(
                    {
                        "peak_index": peak_index,
                        "peak_kind": "fragment",
                        "ion_kind": record.ion_kinds[peak_index],
                        "ion_position": int(record.ion_positions[peak_index]),
                        "mz": _format_value(float(mz)),
                        "theoretical_intensity": _format_value(float(record.theoretical_intensity[peak_index])),
                        "experimental_intensity": _format_value(float(record.experimental_intensity[peak_index])),
                    }
                )
                writer.writerow(row)
                row_count += 1
    return row_count


def _count_summary(records, counter):
    if not records:
        return None
    counts = np.asarray([counter(record) for record in records], dtype=float)
    return {
        "min": float(np.min(counts)),
        "max": float(np.max(counts)),
        "median": float(np.median(counts)),
        "average": float(np.mean(counts)),
    }


def _plot_records(records, input_path, output_path, args):
    plotted = records[: min(MAX_PLOT_RECORDS, len(records))]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(plotted),
        3,
        figsize=(30, 5.2 * len(plotted)),
        squeeze=False,
        gridspec_kw={"width_ratios": [1.0, 2.2, 2.2]},
    )
    fig.subplots_adjust(left=0.055, right=0.99, top=0.92, bottom=0.06, hspace=0.45, wspace=0.16)

    for idx, record in enumerate(plotted):
        _plot_precursor(axes[idx][0], record)
        _plot_fragments(axes[idx][1], record, args.mz_min, args.mz_max)
        _plot_cached_prediction(axes[idx][2], record, args.mz_min, args.mz_max)
        rt = _safe_float(record.retention)
        rt_text = f"{rt:.4f}" if rt is not None else (record.retention or "NA")
        abundance_text = _format_value(record.target_sip_abundance_pct)
        source_text = record.source_file
        if record.record_kind or abundance_text:
            source_text += f" | {record.record_kind or 'record'}"
            if abundance_text:
                source_text += f" {record.sip_label}={abundance_text}%"
        title = (
            f"{_truncate_text(source_text, 92)}\n"
            f"{_truncate_text(record.psm_id)} | RT={rt_text} | z={record.charge} | "
            f"cached iRT={record.predicted_rt:.4f} | "
            f"matched envelopes={matched_envelope_count(record)}/{retained_envelope_count(record)} | "
            f"supported modeled peaks={matched_peak_count(record)}\n"
            f"Peptide={_truncate_text(record.peptide, 72)} | "
            f"Proteins={_truncate_text(record.proteins or 'NA', 86)}"
        )
        axes[idx][0].text(0.0, 1.1, title, transform=axes[idx][0].transAxes, fontsize=10, va="bottom")

    fig.suptitle(
        f"SFI, cached predictions, and real MS/MS from {input_path.name} "
        f"(selected={len(records)}, plotted={len(plotted)}, rank_by={args.rank_by})",
        fontsize=16,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    if args.show:
        plt.show()
    else:
        plt.close(fig)
    return len(plotted)


def main():
    args = parse_args()
    input_path = Path(args.sfi).resolve()
    spectrum_cache_path = Path(args.spectrum_cache).resolve()
    rt_cache_path = Path(args.rt_cache).resolve()
    sample_hdf5_path = Path(args.sample_hdf5).resolve()
    output_path = Path(args.output).resolve()
    tsv_output_path = Path(args.tsv_output).resolve() if args.tsv_output else output_path.with_suffix(".tsv")

    if args.top_k <= 0:
        raise SystemExit("--top-k must be positive.")
    if args.sample_size < 0:
        raise SystemExit("--sample-size must be non-negative.")
    if args.min_matched_envelopes < 0:
        raise SystemExit("--min-matched-envelopes must be non-negative.")
    if args.mz_min is not None and args.mz_max is not None and args.mz_min >= args.mz_max:
        raise SystemExit("--mz-min must be smaller than --mz-max.")
    for label, path in (
        ("SFI", input_path),
        ("spectrum cache", spectrum_cache_path),
        ("RT cache", rt_cache_path),
        ("real-sample HDF5", sample_hdf5_path),
    ):
        if not path.exists():
            raise SystemExit(f"Required {label} input not found: {path}")

    try:
        records, total_sfi_records, sfi_file_count = load_sfi_candidates(
            input_path, args
        )
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(str(error)) from error
    if not records:
        raise SystemExit(f"No SFI spectra records found in {input_path}")

    selected = _select_records(records, args)
    if not selected:
        raise SystemExit("No records matched the requested filters.")
    try:
        attach_predictions(selected, spectrum_cache_path, rt_cache_path)
        attach_raw_spectra(selected, sample_hdf5_path)
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(str(error)) from error

    plotted_count = _plot_records(selected, input_path, output_path, args)
    tsv_row_count = _write_selected_tsv(selected, tsv_output_path)
    matched_envelope_summary = _count_summary(records, matched_envelope_count)
    retained_envelope_summary = _count_summary(records, retained_envelope_count)
    supported_peak_summary = _count_summary(records, matched_peak_count)

    print(f"input_path={input_path}")
    print(f"sfi_files={sfi_file_count}")
    print(f"total_sfi_records={total_sfi_records}")
    print(f"sampled_sfi_records={len(records)}")
    print(f"spectrum_cache={spectrum_cache_path}")
    print(f"rt_cache={rt_cache_path}")
    print(f"real_sample_hdf5={sample_hdf5_path}")
    for label, summary in (
        ("matched_envelope_count", matched_envelope_summary),
        ("retained_envelope_count", retained_envelope_summary),
        ("experimentally_supported_modeled_peak_count", supported_peak_summary),
    ):
        if summary is not None:
            print(
                f"{label}_"
                f"min={_format_value(summary['min'])}\t"
                f"max={_format_value(summary['max'])}\t"
                f"median={_format_value(summary['median'])}\t"
                f"average={_format_value(summary['average'])}"
            )
    print(f"selected_records={len(selected)}")
    print(f"plotted_records={plotted_count}")
    print(f"tsv_peak_rows={tsv_row_count}")
    for record in selected[:plotted_count]:
        print(
            f"plotted_psm={record.psm_id}\t"
            f"RT={record.retention or 'NA'}\t"
            f"PredictedRT={_format_value(record.predicted_rt)}\t"
            f"Charge={record.charge}\t"
            f"Proteins={record.proteins or 'NA'}\t"
            f"matched_envelopes={matched_envelope_count(record)}\t"
            f"retained_envelopes={retained_envelope_count(record)}\t"
            f"supported_modeled_peaks={matched_peak_count(record)}"
        )
    print(f"output_plot={output_path}")
    print(f"output_tsv={tsv_output_path}")


if __name__ == "__main__":
    main()
