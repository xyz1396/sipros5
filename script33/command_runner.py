import os
import subprocess
from logging import Logger


def run_logged_command(cmd: str, logger: Logger, env: dict[str, str] | None = None,
                       env_updates: dict[str, str] | None = None,
                       cwd: str | None = None,
                       cpu_cores: int | None = None) -> None:
    if cpu_cores is None:
        logger.info(f"Running command: {cmd}")
    else:
        unit = "core" if cpu_cores == 1 else "cores"
        logger.info(f"Running process (allocated {cpu_cores} CPU {unit}): {cmd}")
    run_env = os.environ.copy()
    if env is not None:
        run_env.update(env)
    if env_updates is not None:
        for key, value in env_updates.items():
            run_env[key] = str(value)
            logger.debug(f"Process environment: {key}={run_env[key]}")
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, env=run_env, cwd=cwd
        )
        logger.info(output.decode())
    except subprocess.CalledProcessError as exc:
        logger.error(f"Command execution failed: {exc.output.decode()}")
        raise SystemExit(1) from exc
