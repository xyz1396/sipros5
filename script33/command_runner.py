import os
import shlex
import subprocess
import time
from collections.abc import Sequence
from logging import Logger


def quote_argument(value: object) -> str:
    """Quote one argument for the platform shell used below."""
    text = str(value)
    if os.name == "nt":
        return subprocess.list2cmdline([text])
    return shlex.quote(text)


def join_command(arguments) -> str:
    return " ".join(quote_argument(argument) for argument in arguments)


def run_logged_command(cmd: str | Sequence[str], logger: Logger,
                       env: dict[str, str] | None = None,
                       env_updates: dict[str, str] | None = None,
                       cwd: str | None = None,
                       cpu_threads: int | None = None) -> None:
    display_command = cmd if isinstance(cmd, str) else join_command(cmd)
    if cpu_threads is None:
        logger.info(f"Running command: {display_command}")
    else:
        unit = "thread" if cpu_threads == 1 else "threads"
        logger.info(
            f"Running process ({cpu_threads} CPU {unit}): {display_command}"
        )
    run_env = os.environ.copy()
    if env is not None:
        run_env.update(env)
    if env_updates is not None:
        for key, value in env_updates.items():
            run_env[key] = str(value)
            logger.debug(f"Process environment: {key}={run_env[key]}")
    started = time.perf_counter()
    try:
        # Argument arrays bypass cmd.exe on Windows. Besides being safer for
        # quoting, this avoids cmd.exe's 8191-character command-line limit for
        # cross-sample Aerith runs containing many input files.
        use_shell = isinstance(cmd, str)
        output = subprocess.check_output(
            cmd, shell=use_shell, stderr=subprocess.STDOUT, env=run_env, cwd=cwd,
            text=True, encoding="utf-8", errors="replace",
        )
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        decoded = output.rstrip()
        if decoded:
            logger.info(decoded)
        logger.info(
            f"Process completed in {time.perf_counter() - started:.3f} s"
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            f"Process failed after {time.perf_counter() - started:.3f} s: "
            f"{exc.output.decode('utf-8', errors='replace') if isinstance(exc.output, bytes) else exc.output}"
        )
        raise SystemExit(1) from exc
