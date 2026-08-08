from __future__ import annotations

import importlib.util
from pathlib import Path
import struct
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("visualize_experimental_spectra.py")
SPEC = importlib.util.spec_from_file_location("visualize_experimental_spectra", SCRIPT)
VISUALIZE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = VISUALIZE
SPEC.loader.exec_module(VISUALIZE)


class VisualizeSfiTests(unittest.TestCase):
    def test_prediction_cache_uses_normalized_peptide_charge_key(self) -> None:
        self.assertEqual(
            VISUALIZE._peptide_body("R[QMDVVEQMMPGLK]D"),
            "QMDVVEQMMPGLK",
        )

    def test_sfi_psm_id_resolves_real_sample_scan(self) -> None:
        self.assertEqual(
            VISUALIZE._sample_name_and_scan(
                "PanC_20260708_03_DDA_01.12988.1_000.000Pct"
            ),
            ("PanC_20260708_03_DDA_01", 12988),
        )
        self.assertEqual(
            VISUALIZE._sample_name_and_scan(
                "DECOY_PanC_20260708_03_DDA_01.12988.1_050.000Pct"
            ),
            ("PanC_20260708_03_DDA_01", 12988),
        )

    def test_binary_prediction_caches_have_no_missing_key_fallback(self) -> None:
        spectrum_key = b"PEPTIDEK\x1f2"
        with tempfile.TemporaryDirectory() as directory:
            cache_path = Path(directory, "predictions.bin")
            with cache_path.open("wb") as stream:
                stream.write(struct.pack(
                    "<QQQ", VISUALIZE.PREDICTION_CACHE_MAGIC, 7, 1
                ))
                stream.write(struct.pack("<I", len(spectrum_key)))
                stream.write(spectrum_key)
                stream.write(struct.pack(
                    "<BfI",
                    VISUALIZE.PREDICTION_HAS_SPECTRUM |
                    VISUALIZE.PREDICTION_HAS_RT,
                    12.5,
                    1,
                ))
                stream.write(struct.pack("<ffcIi", 500.25, 0.75, b"y", 4, 1))

            spectra, retention = VISUALIZE.load_prediction_cache(
                cache_path, {spectrum_key}
            )
            spectrum = spectra[spectrum_key]
            self.assertAlmostEqual(float(spectrum.mz[0]), 500.25)
            self.assertEqual(spectrum.ion_kinds, ["y"])
            self.assertAlmostEqual(retention[spectrum_key], 12.5)
            with self.assertRaisesRegex(ValueError, "fallback is disabled"):
                VISUALIZE.load_prediction_cache(
                    cache_path, {b"MISSING\x1f2"}
                )


if __name__ == "__main__":
    unittest.main()
