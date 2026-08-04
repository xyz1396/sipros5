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
    def test_rt_cache_key_matches_native_diann_tokenizer(self) -> None:
        self.assertEqual(
            VISUALIZE._rt_token_key("R[QMDVVEQMMPGLK]D"),
            bytes([1, 16, 11, 19, 5, 5, 17, 16, 11, 11, 8, 3, 7, 20, 2]),
        )
        self.assertEqual(
            VISUALIZE._rt_token_key("R[QM~DVVEQMMPGLK]D")[2], 26
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
        rt_key = bytes([1, 8, 17, 8, 14, 6, 19, 17, 20, 2])
        with tempfile.TemporaryDirectory() as directory:
            spectrum_path = Path(directory, "predictions.spectrum")
            rt_path = Path(directory, "predictions.rt")
            with spectrum_path.open("wb") as stream:
                stream.write(struct.pack(
                    "<QQQ", VISUALIZE.SPECTRUM_CACHE_MAGIC, 7, 1
                ))
                stream.write(struct.pack("<I", len(spectrum_key)))
                stream.write(spectrum_key)
                stream.write(struct.pack("<I", 1))
                stream.write(struct.pack("<ffcIi", 500.25, 0.75, b"y", 4, 1))
            with rt_path.open("wb") as stream:
                stream.write(struct.pack("<QQQ", VISUALIZE.RT_CACHE_MAGIC, 9, 1))
                stream.write(struct.pack("<I", len(rt_key)))
                stream.write(rt_key)
                stream.write(struct.pack("<f", 12.5))

            spectrum = VISUALIZE.load_spectrum_cache(
                spectrum_path, {spectrum_key}
            )[spectrum_key]
            self.assertAlmostEqual(float(spectrum.mz[0]), 500.25)
            self.assertEqual(spectrum.ion_kinds, ["y"])
            self.assertAlmostEqual(
                VISUALIZE.load_rt_cache(rt_path, {rt_key})[rt_key], 12.5
            )
            with self.assertRaisesRegex(ValueError, "fallback is disabled"):
                VISUALIZE.load_spectrum_cache(
                    spectrum_path, {b"MISSING\x1f2"}
                )
            with self.assertRaisesRegex(ValueError, "fallback is disabled"):
                VISUALIZE.load_rt_cache(rt_path, {b"missing"})


if __name__ == "__main__":
    unittest.main()
