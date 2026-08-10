import os
import subprocess
import time
from logging import Logger


def run_logged_command(cmd: str, logger: Logger, env: dict[str, str] | None = None,
                       env_updates: dict[str, str] | None = None,
                       cwd: str | None = None,
                       cpu_threads: int | None = None) -> None:
    if cpu_threads is None:
        logger.info(f"Running command: {cmd}")
    else:
        unit = "thread" if cpu_threads == 1 else "threads"
        logger.info(
            f"Running process ({cpu_threads} CPU {unit}): {cmd}"
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
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, env=run_env, cwd=cwd
        )
        decoded = output.decode().rstrip()
        if decoded:
            logger.info(decoded)
        logger.info(
            f"Process completed in {time.perf_counter() - started:.3f} s"
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            f"Process failed after {time.perf_counter() - started:.3f} s: "
            f"{exc.output.decode()}"
        )
        raise SystemExit(1) from exc
