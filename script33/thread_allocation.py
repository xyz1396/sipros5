"""Helpers for keeping nested workflow parallelism within one CPU budget."""

from __future__ import annotations

from dataclasses import dataclass
import os


MIN_SIPROS_THREADS = 8


@dataclass(frozen=True)
class ThreadAllocation:
    """A concurrency limit and the thread budget assigned to each task."""

    worker_count: int
    task_threads: tuple[int, ...]

    @property
    def peak_threads(self) -> int:
        """Maximum threads used when the first wave of tasks is running."""
        return sum(self.task_threads[:self.worker_count])


def available_cpu_count() -> int:
    """Return CPUs available to this process, respecting CPU affinity."""
    process_cpu_count = getattr(os, "process_cpu_count", None)
    if process_cpu_count is not None:
        count = process_cpu_count()
        if count:
            return max(1, count)
    try:
        affinity_count = len(os.sched_getaffinity(0))
        if affinity_count:
            return affinity_count
    except (AttributeError, OSError, NotImplementedError):
        pass
    return max(1, os.cpu_count() or 1)


def effective_thread_count(requested: int, available: int | None = None) -> int:
    """Resolve ``0`` to all available CPUs and cap explicit over-allocation."""
    if requested < 0:
        raise ValueError("Thread count must be non-negative")
    cpu_count = available_cpu_count() if available is None else max(1, available)
    if requested == 0:
        return cpu_count
    return min(requested, cpu_count)


def allocate_threads(total_threads: int, task_count: int, *,
                     minimum_threads_per_task: int = 1) -> ThreadAllocation:
    """Split a total thread budget across concurrently runnable tasks.

    The first wave receives an even split of the total budget. Queued tasks use
    the smaller base allocation, so replacing a completed first-wave task can
    never increase live concurrency above ``total_threads``.

    ``minimum_threads_per_task`` limits the number of simultaneous workers. If
    the total budget is smaller than that minimum, one task receives the whole
    budget rather than oversubscribing beyond the user's limit.
    """
    if total_threads <= 0:
        raise ValueError("total_threads must be positive")
    if task_count < 0:
        raise ValueError("task_count must be non-negative")
    if minimum_threads_per_task <= 0:
        raise ValueError("minimum_threads_per_task must be positive")
    if task_count == 0:
        return ThreadAllocation(worker_count=0, task_threads=())

    effective_minimum = min(total_threads, minimum_threads_per_task)
    worker_count = min(
        task_count, max(1, total_threads // effective_minimum)
    )
    threads_per_task, remainder = divmod(total_threads, worker_count)
    first_wave = tuple(
        threads_per_task + (index < remainder)
        for index in range(worker_count)
    )
    queued = (threads_per_task,) * (task_count - worker_count)
    task_threads = first_wave + queued
    return ThreadAllocation(worker_count=worker_count, task_threads=task_threads)


def thread_env_updates(thread_count: int) -> dict[str, str]:
    """Environment limits understood by the native, Go, and .NET tools."""
    if thread_count <= 0:
        raise ValueError("thread_count must be positive")
    value = str(thread_count)
    return {
        "OMP_NUM_THREADS": value,
        "OMP_THREAD_LIMIT": value,
        "OMP_MAX_ACTIVE_LEVELS": "1",
        "OMP_DYNAMIC": "FALSE",
        # A task gets one shared quota, not a separate quota for every nested
        # runtime.  Sipros uses OpenMP, Philosopher uses Go, and Raxport uses
        # .NET; limiting incidental BLAS pools to one prevents q-by-q nesting.
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "GOMAXPROCS": value,
        "DOTNET_PROCESSOR_COUNT": value,
        "COMPlus_ProcessorCount": value,
    }


def configure_process_worker(thread_count: int = 1) -> None:
    """Limit native libraries inside a Python process-pool worker."""
    os.environ.update(thread_env_updates(thread_count))
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=thread_count)
    except ImportError:
        pass
