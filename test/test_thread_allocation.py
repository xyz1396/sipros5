#!/usr/bin/env python3
"""Regression tests for the workflow-wide CPU-thread budget."""

from __future__ import annotations

import logging
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
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
        return workflow

    def test_workflows_write_config_without_configs_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            templates = root / "templates"
            templates.mkdir()
            (templates / "Regular.cfg").write_text(
                "Search_Name = old\nSIP_Element = C\nSIP_Element_Isotope = 13\n"
            )
            (templates / "SIP.cfg").write_text(
                "Search_Name = old\nSIP_Element = C\nSIP_Element_Isotope = 13\n"
            )
            for element, filename in (("R", "Regular.cfg"), ("N15", "SIP.cfg")):
                with self.subTest(element=element):
                    output = root / element
                    workflow = self.make_search(str(output))
                    workflow.element = element
                    workflow.configTemplatePath = str(templates)
                    workflow.fastaPath = "target.fasta"

                    config = workflow.write_workflow_config()

                    self.assertEqual(config, str(output / filename))
                    self.assertTrue(Path(config).is_file())
                    self.assertFalse((output / "configs").exists())

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

            workflow.search_spectra_samples("config.cfg", "spectra")

            self.assertEqual(sorted(threads for _, threads in captured), [8, 8, 8])
            for command, threads in captured:
                match = re.search(r"(?:^| )-t (\d+)(?: |$)", command)
                self.assertIsNotNone(match)
                self.assertEqual(int(match.group(1)), threads)

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

            workflow.search_spectra_samples("config.cfg", "spectra")

            self.assertEqual(sorted(captured), [8, 8, 8])

    def test_fasta_target_and_decoy_share_one_budget(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["sample"]
            workflow.hdf5_paths = {"sample": f"{output}/sample.h5"}
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            workflow.direct_sip_args = lambda: ""
            captured: list[int] = []
            workflow.run_command_sipros = (
                lambda _command, threads: captured.append(threads)
            )
            workflow.merge_pin_files = lambda *_arguments: None

            workflow.sipros_search("config.cfg")

            self.assertEqual(sorted(captured), [8, 8])
            self.assertEqual(sum(captured), 16)

    def test_fasta_search_overwrites_existing_pin(self) -> None:
        with tempfile.TemporaryDirectory() as output:
            workflow = self.make_search(output)
            workflow.base_names = ["sample"]
            workflow.hdf5_paths = {"sample": f"{output}/sample.h5"}
            workflow.element = "R"
            workflow.fastaPath = "target.fasta"
            workflow.decoyPath = "decoy.fasta"
            workflow.direct_sip_args = lambda: ""
            completed = Path(output) / "sample" / "sample_target.pin"
            completed.parent.mkdir(parents=True)
            completed.write_bytes(b"x" * (500 * 1024 + 1))
            captured: list[int] = []
            merged: list[tuple[str, str, str]] = []
            workflow.run_command_sipros = (
                lambda _command, threads: captured.append(threads)
            )
            workflow.merge_pin_files = (
                lambda *arguments: merged.append(arguments)
            )

            workflow.sipros_search("config.cfg")

            self.assertEqual(captured, [8, 8])
            self.assertEqual(len(merged), 1)

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


if __name__ == "__main__":
    unittest.main()
