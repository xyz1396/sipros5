#!/usr/bin/env python3
"""Regression tests for the workflow-wide CPU-thread budget."""

from __future__ import annotations

import logging
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "script33"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import assembly as assembly_module
import command_runner as command_runner_module
import filter as filter_module
import search as search_module
import thread_allocation as allocation_module


def null_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    logger.propagate = False
    return logger


class CommandLoggingTests(unittest.TestCase):
    def test_process_log_reports_cores_without_environment_spam(self) -> None:
        logger = mock.Mock()
        with mock.patch.object(
                command_runner_module.subprocess,
                "check_output",
                return_value=b"complete\n",
        ) as check_output:
            command_runner_module.run_logged_command(
                "tool --input sample",
                logger,
                env_updates={"OMP_NUM_THREADS": "6", "GOMAXPROCS": "6"},
                cpu_cores=6,
            )

        self.assertEqual(
            logger.info.call_args_list[0].args[0],
            "Running process (allocated 6 CPU cores): tool --input sample",
        )
        info_messages = [call.args[0] for call in logger.info.call_args_list]
        self.assertFalse(any(message.startswith("Set ") for message in info_messages))
        process_env = check_output.call_args.kwargs["env"]
        self.assertEqual(process_env["OMP_NUM_THREADS"], "6")
        self.assertEqual(process_env["GOMAXPROCS"], "6")


class ThreadAllocationTests(unittest.TestCase):
    def test_exact_budget_splits(self) -> None:
        cases = {
            (8, 0): (0, ()),
            (8, 1): (1, (8,)),
            (8, 3): (3, (3, 3, 2)),
            (8, 10): (8, (1,) * 10),
            (16, 3): (3, (6, 5, 5)),
            (16, 10): (10, (2,) * 6 + (1,) * 4),
            (16, 32): (16, (1,) * 32),
            (32, 3): (3, (11, 11, 10)),
            (32, 10): (10, (4, 4) + (3,) * 8),
        }
        for arguments, expected in cases.items():
            with self.subTest(total_threads=arguments[0], task_count=arguments[1]):
                allocation = allocation_module.allocate_threads(*arguments)
                self.assertEqual(
                    (allocation.worker_count, allocation.task_threads), expected
                )

    def test_allocation_invariants(self) -> None:
        for total_threads in range(1, 65):
            for task_count in range(0, 100):
                allocation = allocation_module.allocate_threads(
                    total_threads, task_count
                )
                self.assertEqual(len(allocation.task_threads), task_count)
                self.assertEqual(
                    allocation.worker_count, min(total_threads, task_count)
                )
                if task_count == 0:
                    self.assertEqual(allocation.peak_threads, 0)
                    continue
                self.assertTrue(all(value >= 1 for value in allocation.task_threads))
                self.assertEqual(allocation.peak_threads, total_threads)
                if task_count <= total_threads:
                    self.assertEqual(sum(allocation.task_threads), total_threads)
                    self.assertLessEqual(
                        max(allocation.task_threads) - min(allocation.task_threads), 1
                    )
                else:
                    self.assertEqual(set(allocation.task_threads), {1})

    def test_invalid_allocations(self) -> None:
        for total_threads, task_count in ((0, 1), (-1, 1), (1, -1)):
            with self.subTest(total_threads=total_threads, task_count=task_count):
                with self.assertRaises(ValueError):
                    allocation_module.allocate_threads(total_threads, task_count)
        with self.assertRaises(ValueError):
            allocation_module.allocate_threads(
                8, 1, minimum_threads_per_task=0
            )

    def test_minimum_eight_thread_process_allocations(self) -> None:
        cases = {
            (8, 10): (1, (8,) * 10),
            (15, 10): (1, (15,) * 10),
            (16, 10): (2, (8,) * 10),
            (17, 10): (2, (9,) + (8,) * 9),
            (25, 10): (3, (9,) + (8,) * 9),
            (32, 10): (4, (8,) * 10),
        }
        for arguments, expected in cases.items():
            with self.subTest(total_threads=arguments[0], task_count=arguments[1]):
                allocation = allocation_module.allocate_threads(
                    *arguments, minimum_threads_per_task=8
                )
                self.assertEqual(
                    (allocation.worker_count, allocation.task_threads), expected
                )
                self.assertEqual(allocation.peak_threads, arguments[0])
                self.assertGreaterEqual(min(allocation.task_threads), 8)

        for total_threads in range(1, 8):
            below_minimum = allocation_module.allocate_threads(
                total_threads, 10, minimum_threads_per_task=8
            )
            self.assertEqual(below_minimum.worker_count, 1)
            self.assertEqual(
                below_minimum.task_threads, (total_threads,) * 10
            )

    def test_minimum_thread_allocation_invariants(self) -> None:
        for total_threads in range(1, 65):
            for task_count in range(0, 100):
                allocation = allocation_module.allocate_threads(
                    total_threads,
                    task_count,
                    minimum_threads_per_task=8,
                )
                if task_count == 0:
                    self.assertEqual(allocation.worker_count, 0)
                    continue
                expected_minimum = min(total_threads, 8)
                self.assertEqual(
                    allocation.worker_count,
                    min(
                        task_count,
                        max(1, total_threads // expected_minimum),
                    ),
                )
                self.assertGreaterEqual(
                    min(allocation.task_threads), expected_minimum
                )
                largest_live_quotas = sorted(
                    allocation.task_threads, reverse=True
                )[:allocation.worker_count]
                self.assertEqual(sum(largest_live_quotas), total_threads)

    def test_available_cpu_count_prefers_process_count(self) -> None:
        with mock.patch.object(
                allocation_module.os, "process_cpu_count", return_value=6,
                create=True), mock.patch.object(
                    allocation_module.os, "sched_getaffinity", return_value=set(range(12))
                ):
            self.assertEqual(allocation_module.available_cpu_count(), 6)

    def test_available_cpu_count_falls_back_to_affinity_and_cpu_count(self) -> None:
        with mock.patch.object(
                allocation_module.os, "process_cpu_count", return_value=None,
                create=True), mock.patch.object(
                    allocation_module.os, "sched_getaffinity", return_value=set(range(7))
                ):
            self.assertEqual(allocation_module.available_cpu_count(), 7)

        with mock.patch.object(
                allocation_module.os, "process_cpu_count", return_value=None,
                create=True), mock.patch.object(
                    allocation_module.os, "sched_getaffinity", side_effect=OSError
                ), mock.patch.object(allocation_module.os, "cpu_count", return_value=5):
            self.assertEqual(allocation_module.available_cpu_count(), 5)

        with mock.patch.object(
                allocation_module.os, "process_cpu_count", return_value=None,
                create=True), mock.patch.object(
                    allocation_module.os, "sched_getaffinity",
                    side_effect=NotImplementedError
                ), mock.patch.object(allocation_module.os, "cpu_count", return_value=None):
            self.assertEqual(allocation_module.available_cpu_count(), 1)

    def test_effective_thread_count_caps_to_available_cpus(self) -> None:
        self.assertEqual(allocation_module.effective_thread_count(0, 64), 64)
        self.assertEqual(allocation_module.effective_thread_count(0, 4), 4)
        self.assertEqual(allocation_module.effective_thread_count(8, 64), 8)
        self.assertEqual(allocation_module.effective_thread_count(16, 64), 16)
        self.assertEqual(allocation_module.effective_thread_count(32, 64), 32)
        self.assertEqual(allocation_module.effective_thread_count(32, 16), 16)
        with self.assertRaises(ValueError):
            allocation_module.effective_thread_count(-1, 16)

    def test_child_environment_caps_each_runtime(self) -> None:
        environment = allocation_module.thread_env_updates(6)
        self.assertEqual(environment["OMP_NUM_THREADS"], "6")
        self.assertEqual(environment["OMP_THREAD_LIMIT"], "6")
        self.assertEqual(environment["GOMAXPROCS"], "6")
        self.assertEqual(environment["DOTNET_PROCESSOR_COUNT"], "6")
        self.assertEqual(environment["COMPlus_ProcessorCount"], "6")
        self.assertEqual(environment["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual(environment["OMP_MAX_ACTIVE_LEVELS"], "1")
        self.assertEqual(environment["OMP_DYNAMIC"], "FALSE")

    @unittest.skipUnless(Path("/proc/self/status").exists(), "Linux /proc required")
    def test_pandas_loads_native_pool_after_one_thread_limit(self) -> None:
        environment = os.environ.copy()
        for variable in allocation_module.thread_env_updates(1):
            environment.pop(variable, None)
        code = (
            "import sys; "
            f"sys.path.insert(0, {str(SCRIPT_DIR)!r}); "
            "import filter; "
            "import numpy as np; "
            "values = np.ones((256, 256)); values @ values; "
            "print(next(line.split()[1] for line in open('/proc/self/status') "
            "if line.startswith('Threads:')))"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        self.assertLessEqual(int(result.stdout.strip()), 2, result.stderr)


class WorkflowAllocationTests(unittest.TestCase):
    def make_search(self, output: str, threads: int = 16):
        workflow = object.__new__(search_module.search)
        workflow.threadNumber = threads
        workflow.logger = null_logger(f"search-allocation-{id(workflow)}")
        workflow.outPutPath = output
        workflow.siprosPath = "sipros"
        workflow.toleranceMS1 = 0.01
        workflow.toleranceMS2 = 0.02
        workflow.topPsmsPerScan = 8
        workflow.ptms = None
        workflow.fixedPtms = None
        workflow.maxPtmCount = None
        workflow.decoyPrefix = "Decoy_"
        return workflow

    def test_regular_decoy_preserves_protein_n_terminal_residue(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            fasta = Path(output, "target.fasta")
            fasta.write_text(">met\nMABCDEK\n>nonmet\nABCDEK\n")
            workflow = self.make_search(output)
            workflow.element = "R"
            workflow.fastaPath = str(fasta)
            workflow.decoyPath = str(Path(output, "decoy.fasta"))

            workflow.reverse_fasta_sequences()

            self.assertEqual(
                Path(workflow.decoyPath).read_text(),
                ">Decoy_met\nMKEDCBA\n>Decoy_nonmet\nAKEDCB\n",
            )

    def test_sip_decoy_keeps_legacy_full_reversal(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            fasta = Path(output, "target.fasta")
            fasta.write_text(">protein\nMABCDEK\n")
            workflow = self.make_search(output)
            workflow.element = "C13"
            workflow.fastaPath = str(fasta)
            workflow.decoyPath = str(Path(output, "decoy.fasta"))

            workflow.reverse_fasta_sequences()

            self.assertEqual(
                Path(workflow.decoyPath).read_text(),
                ">Decoy_protein\nKEDCBAM\n",
            )

    def test_spectra_search_divides_cli_threads_between_samples(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["one", "two", "three"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            captured: list[tuple[str, int]] = []
            workflow.run_command_sipros = (
                lambda command, threads: captured.append((command, threads))
            )

            workflow.search_spectra_samples("spectra")

            self.assertEqual(sorted(threads for _, threads in captured), [8, 8, 8])
            for command, threads in captured:
                match = re.search(r"(?:^| )-t (\d+)(?: |$)", command)
                self.assertIsNotNone(match)
                self.assertEqual(int(match.group(1)), threads)
                self.assertNotIn(" -c ", command)

    def test_spectra_search_overwrites_existing_pin(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["done", "two", "three"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            completed = Path(output) / "done" / "done.pin"
            completed.parent.mkdir(parents=True)
            completed.write_bytes(b"x" * (500 * 1024 + 1))
            captured: list[int] = []
            workflow.run_command_sipros = (
                lambda _command, threads: captured.append(threads)
            )

            workflow.search_spectra_samples("spectra")

            self.assertEqual(sorted(captured), [8, 8, 8])

    def test_fasta_target_and_decoy_share_one_budget(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output, threads=24)
            workflow.base_names = ["one", "two", "three"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            captured: list[tuple[str, int]] = []
            workflow.run_command_sipros = (
                lambda command, threads: captured.append((command, threads))
            )
            merged: list[tuple[str, str, str]] = []
            workflow.merge_pin_files = lambda *arguments: merged.append(arguments)

            workflow.sipros_search()

            self.assertEqual(len(captured), 8)
            prepare_calls = captured[:2]
            search_calls = captured[2:]
            self.assertEqual(
                sorted(threads for _, threads in prepare_calls), [12, 12]
            )
            self.assertEqual(
                sum(threads for _, threads in prepare_calls), 24
            )
            self.assertEqual(
                sorted(threads for _, threads in search_calls), [12] * 6
            )
            cache_paths: set[str] = set()
            output_paths: set[str] = set()
            for command, _ in prepare_calls:
                arguments = shlex.split(command)
                self.assertIn("--prepare-only", arguments)
                self.assertNotIn("-f", arguments)
                self.assertNotIn("-o", arguments)
                cache_index = arguments.index("--fragment-index-cache")
                cache_paths.add(arguments[cache_index + 1])
                self.assertEqual(
                    Path(arguments[cache_index + 1]).parent, Path(output)
                )
            for command, _ in search_calls:
                self.assertNotIn(" -c ", command)
                self.assertIn("--tolerance-ms1 0.01", command)
                self.assertIn("--tolerance-ms2 0.02", command)
                arguments = shlex.split(command)
                self.assertNotIn("--prepare-only", arguments)
                scan_files = [
                    arguments[index + 1]
                    for index, argument in enumerate(arguments)
                    if argument == "-f"
                ]
                self.assertEqual(len(scan_files), 1)
                self.assertIn(scan_files[0], workflow.hdf5_paths.values())
                self.assertIn("--pin-output", arguments)
                pin_index = arguments.index("--pin-output")
                self.assertIn(
                    arguments[pin_index + 1],
                    {f"{name}_target.pin" for name in workflow.base_names}
                    | {f"{name}_decoy.pin" for name in workflow.base_names},
                )
                self.assertIn("--precursor-source", arguments)
                source_index = arguments.index("--precursor-source")
                self.assertEqual(arguments[source_index + 1], "ms1-neighborhood")
                cache_index = arguments.index("--fragment-index-cache")
                cache_paths.add(arguments[cache_index + 1])
                output_index = arguments.index("-o")
                output_paths.add(arguments[output_index + 1])
                self.assertNotIn("--ptm", arguments)
                self.assertNotIn("--fixed-ptm", arguments)
                self.assertNotIn("--max-ptm-count", arguments)
            self.assertEqual(
                {Path(path).name for path in cache_paths},
                {"target.sfi", "decoy.sfi"},
            )
            self.assertEqual(
                output_paths,
                {str(Path(output) / name) for name in workflow.base_names},
            )
            self.assertEqual(len(merged), 3)
            self.assertFalse((Path(output) / "regular_fasta_search").exists())

    def test_regular_fasta_pairs_split_odd_budget_for_every_sample(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output, threads=5)
            workflow.base_names = ["one", "two"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            captured: list[tuple[str, int]] = []
            workflow.run_command_sipros = (
                lambda command, threads: captured.append((command, threads))
            )
            merged: list[tuple[str, str, str]] = []
            workflow.merge_pin_files = lambda *arguments: merged.append(arguments)

            workflow.sipros_search()

            search_calls = captured[2:]
            self.assertEqual(len(search_calls), 4)
            per_sample: dict[str, dict[str, int]] = {
                name: {} for name in workflow.base_names
            }
            for command, threads in search_calls:
                arguments = shlex.split(command)
                scan = arguments[arguments.index("-f") + 1]
                sample = Path(scan).stem
                label = arguments[arguments.index("--pin-label") + 1]
                per_sample[sample][label] = threads
            self.assertEqual(
                per_sample,
                {
                    "one": {"1": 3, "-1": 2},
                    "two": {"1": 3, "-1": 2},
                },
            )
            self.assertEqual(len(merged), 2)

    def test_regular_fasta_finishes_each_pair_before_next_sample(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output, threads=8)
            workflow.base_names = ["one", "two"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            barriers = {
                name: threading.Barrier(2) for name in workflow.base_names
            }
            lock = threading.Lock()
            active_counts = {name: 0 for name in workflow.base_names}
            maximum_active_samples = 0

            def fake_run(command: str, _threads: int) -> None:
                nonlocal maximum_active_samples
                arguments = shlex.split(command)
                if "--prepare-only" in arguments:
                    return
                scan = arguments[arguments.index("-f") + 1]
                sample = Path(scan).stem
                with lock:
                    active_counts[sample] += 1
                    maximum_active_samples = max(
                        maximum_active_samples,
                        sum(count > 0 for count in active_counts.values()),
                    )
                barriers[sample].wait(timeout=2)
                with lock:
                    active_counts[sample] -= 1

            workflow.run_command_sipros = fake_run
            workflow.merge_pin_files = lambda *_arguments: None

            workflow.sipros_search()

            self.assertEqual(maximum_active_samples, 1)
            self.assertTrue(all(barrier.broken is False for barrier in barriers.values()))

    def test_fasta_search_passes_custom_ptms_to_target_and_decoy(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["sample"]
            workflow.hdf5_paths = {"sample": f"{output}/sample.h5"}
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            workflow.ptms = ["phosphorylation", "!"]
            workflow.fixedPtms = ["none"]
            workflow.maxPtmCount = 2
            captured: list[str] = []
            workflow.run_command_sipros = (
                lambda command, _threads: captured.append(command)
            )
            workflow.merge_pin_files = lambda *_arguments: None

            workflow.sipros_search()

            self.assertEqual(len(captured), 4)
            for command in captured:
                arguments = shlex.split(command)
                ptms = [
                    arguments[index + 1]
                    for index, argument in enumerate(arguments)
                    if argument == "--ptm"
                ]
                self.assertEqual(ptms, ["phosphorylation", "!"])
                fixed_ptms = [
                    arguments[index + 1]
                    for index, argument in enumerate(arguments)
                    if argument == "--fixed-ptm"
                ]
                self.assertEqual(fixed_ptms, ["none"])
                max_index = arguments.index("--max-ptm-count")
                self.assertEqual(arguments[max_index + 1], "2")

    def test_sip_fasta_keeps_per_sample_legacy_commands(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output, threads=24)
            workflow.base_names = ["one", "two"]
            workflow.hdf5_paths = {
                name: f"{output}/{name}.h5" for name in workflow.base_names
            }
            workflow.element = "C13"
            workflow.sipRange = "0-5"
            workflow.step = "1"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            captured: list[tuple[str, int]] = []
            workflow.run_command_sipros = (
                lambda command, threads: captured.append((command, threads))
            )
            workflow.merge_pin_files = lambda *_arguments: None

            workflow.sipros_search()

            self.assertEqual(len(captured), 4)
            self.assertEqual(sorted(threads for _, threads in captured), [8] * 4)
            for command, _ in captured:
                arguments = shlex.split(command)
                self.assertEqual(arguments.count("-f"), 1)
                self.assertIn("--pin-output", arguments)
                self.assertIn("-a", arguments)
                self.assertIn("-b", arguments)
                self.assertIn("-s", arguments)
                self.assertNotIn("--fragment-index-cache", arguments)
                self.assertNotIn("--precursor-source", arguments)

    def test_generated_spectra_library_receives_fixed_ptms(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.spectraDir = None
            workflow.psmTsv = "input psm.tsv"
            workflow.unlabeledInput = "unlabeled.h5"
            workflow.generatedSpectraDir = f"{output}/spectra"
            workflow.sipRange = "1-5"
            workflow.step = "1"
            workflow.element = "C13"
            workflow.fixedPtms = ["default", "carbamidomethyl"]
            workflow.dryrun = False
            workflow.resolve_or_convert_unlabeled_hdf5 = lambda: "unlabeled.h5"
            captured: list[str] = []

            def fake_run(command: str, threads: int) -> None:
                captured.append(command)
                Path(workflow.generatedSpectraDir, "library.h5").touch()

            workflow.run_command = fake_run

            self.assertEqual(
                workflow.generate_or_reuse_spectra_library(),
                workflow.generatedSpectraDir,
            )
            self.assertEqual(len(captured), 1)
            arguments = shlex.split(captured[0])
            fixed_ptms = [
                arguments[index + 1]
                for index, argument in enumerate(arguments)
                if argument == "--fixed-ptm"
            ]
            self.assertEqual(fixed_ptms, ["default", "carbamidomethyl"])

    def test_reused_spectra_library_rejects_fixed_ptms(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            spectra_dir = Path(output, "spectra")
            spectra_dir.mkdir()
            Path(spectra_dir, "library.h5").touch()
            workflow = self.make_search(output)
            workflow.spectraDir = str(spectra_dir)
            workflow.fixedPtms = ["none"]
            workflow.logger = mock.Mock()

            with self.assertRaises(SystemExit):
                workflow.generate_or_reuse_spectra_library()

            workflow.logger.error.assert_called_once()
            self.assertIn(
                "chemistry metadata is authoritative",
                workflow.logger.error.call_args.args[0],
            )

    def test_fasta_search_overwrites_existing_pin(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["sample"]
            workflow.hdf5_paths = {"sample": f"{output}/sample.h5"}
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            completed = Path(output) / "sample" / "sample_target.pin"
            completed.parent.mkdir(parents=True)
            completed.write_bytes(b"x" * (500 * 1024 + 1))
            captured: list[int] = []
            merged: list[tuple[str, str, str]] = []

            def fake_run(command: str, threads: int) -> None:
                captured.append(threads)
                arguments = shlex.split(command)
                if "--prepare-only" in arguments:
                    return
                output_index = arguments.index("-o")
                pin_index = arguments.index("--pin-output")
                label_index = arguments.index("--pin-label")
                pin_path = (
                    Path(arguments[output_index + 1]) /
                    arguments[pin_index + 1]
                )
                pin_path.write_bytes(
                    b"new target"
                    if arguments[label_index + 1] == "1"
                    else b"new decoy"
                )

            workflow.run_command_sipros = fake_run
            workflow.merge_pin_files = (
                lambda *arguments: merged.append(arguments)
            )

            workflow.sipros_search()

            self.assertEqual(captured, [8, 8, 8, 8])
            self.assertEqual(len(merged), 1)
            self.assertEqual(completed.read_bytes(), b"new target")
            self.assertEqual(
                (completed.parent / "sample_decoy.pin").read_bytes(),
                b"new decoy",
            )
            self.assertFalse((Path(output) / "regular_fasta_search").exists())

    def test_raxport_single_file_uses_dotnet_quota_not_dash_j(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.raxportPath = "raxport"
            workflow.nPrecursor = 6
            expected = Path(output) / "sample.h5"
            captured: list[tuple[str, int]] = []

            def fake_run(command, _environment, threads):
                captured.append((command, threads))
                expected.touch()

            workflow.run_command = fake_run
            workflow.run_command_raxport(
                "sample.raw", output, str(expected), threads=6
            )

            self.assertEqual(captured[0][1], 6)
            self.assertNotIn(" -j ", captured[0][0])

    def test_filter_divides_percolator_threads(self) -> None:
        workflow = object.__new__(filter_module.filter)
        workflow.baseNames = ["one", "two", "three"]
        workflow.outPutPath = "/output"
        workflow.percolatorPath = "percolator"
        workflow.threadNumber = 16
        workflow.logger = null_logger(f"filter-allocation-{id(workflow)}")
        workflow.ignorePCT = False
        workflow.dryrun = False
        captured: list[tuple[str, dict[str, str], int]] = []

        def fake_logged_command(command, _logger, **kwargs):
            captured.append((
                command,
                kwargs["env_updates"],
                kwargs["cpu_cores"],
            ))

        with mock.patch.object(
                filter_module, "run_logged_command", side_effect=fake_logged_command):
            workflow.run()

        command_threads = sorted(
            int(re.search(r"--num-threads (\d+)", command).group(1))
            for command, _, _ in captured
        )
        self.assertEqual(command_threads, [8, 8, 8])
        self.assertEqual(
            sorted(int(environment["OMP_NUM_THREADS"])
                   for _, environment, _ in captured),
            [8, 8, 8],
        )
        self.assertEqual(sorted(cores for _, _, cores in captured), [8, 8, 8])

    def test_philosopher_jobs_receive_balanced_gomaxprocs(self) -> None:
        workflow = object.__new__(assembly_module.assembly)
        workflow.threadNumber = 16
        workflow.logger = null_logger(f"assembly-allocation-{id(workflow)}")
        captured: list[tuple[str, str, int]] = []
        workflow.run_command = (
            lambda command, path, threads: captured.append((command, path, threads))
        )

        allocation = workflow.run_parallel_commands(
            ["one", "two", "three"], ["/one", "/two", "/three"],
            "Philosopher test",
        )

        self.assertEqual(allocation.worker_count, 3)
        self.assertEqual(allocation.peak_threads, 16)
        self.assertEqual(sorted(threads for _, _, threads in captured), [5, 5, 6])


class PepXmlModificationTests(unittest.TestCase):
    def test_protein_n_terminal_acetylation_survives_pepxml_conversion(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            sample_dir = Path(output, "sample")
            sample_dir.mkdir()
            workflow = object.__new__(assembly_module.assembly)
            workflow.outputPath = output
            workflow.decoyPrefix = "Decoy_"
            workflow.logger = null_logger(f"assembly-pepxml-{id(workflow)}")
            psms = assembly_module.pd.DataFrame([
                {
                    "ScanNr": 101,
                    "parentCharges": 2,
                    "SpecId": "sample.101.1",
                    "ExpMass": 1000.0,
                    "retentiontime": 12.5,
                    "Peptide": "[%PEPTIDE]K",
                    "Proteins": "{protein_one}",
                    "massErrors": 0.001,
                    "missCleavageSiteNumbers": 0,
                    "ranks": 1,
                    "posterior_error_prob": 0.01,
                    "score": 2.0,
                },
                {
                    "ScanNr": 102,
                    "parentCharges": 2,
                    "SpecId": "sample.102.1",
                    "ExpMass": 900.0,
                    "retentiontime": 13.5,
                    # '%' after K is a residue acetylation, not the special
                    # protein-N-terminal proteoform.
                    "Peptide": "[AK%]",
                    "Proteins": "{protein_two}",
                    "massErrors": 0.002,
                    "missCleavageSiteNumbers": 0,
                    "ranks": 1,
                    "posterior_error_prob": 0.02,
                    "score": 1.5,
                },
            ])

            workflow.dataframe_to_pepxml(psms, "sample")

            tree = assembly_module.etree.parse(
                str(sample_dir / "sample.pep.xml")
            )
            terminal_modifications = tree.xpath(
                "//search_summary/terminal_modification"
            )
            self.assertEqual(len(terminal_modifications), 1)
            self.assertEqual(
                terminal_modifications[0].get("massdiff"), "42.010565"
            )
            self.assertEqual(terminal_modifications[0].get("terminus"), "N")
            self.assertEqual(
                terminal_modifications[0].get("protein_terminus"), "Y"
            )
            summary_children = [
                child.tag for child in tree.xpath("//search_summary")[0]
            ]
            self.assertLess(
                summary_children.index("terminal_modification"),
                summary_children.index("parameter"),
            )

            acetyl_hit = tree.xpath("//search_hit[@peptide='PEPTIDE']")[0]
            modification_info = acetyl_hit.find("modification_info")
            self.assertIsNotNone(modification_info)
            self.assertEqual(
                modification_info.get("modified_peptide"), "n[43]PEPTIDE"
            )
            self.assertAlmostEqual(
                float(modification_info.get("mod_nterm_mass")),
                43.018390,
                places=6,
            )

            lysine_acetyl_hit = tree.xpath("//search_hit[@peptide='AK']")[0]
            self.assertIsNone(lysine_acetyl_hit.find("modification_info"))


if __name__ == "__main__":
    unittest.main()
