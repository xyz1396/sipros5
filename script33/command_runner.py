import os
import subprocess
from logging import Logger


def run_logged_command(cmd: str, logger: Logger, env: dict[str, str] | None = None,
                       env_updates: dict[str, str] | None = None, cwd: str | None = None) -> None:
    logger.info(f"Running command: {cmd}")
    run_env = os.environ.copy()
    if env is not None:
        run_env.update(env)
    if env_updates is not None:
        for key, value in env_updates.items():
            run_env[key] = str(value)
            logger.info(f"Set {key} to {run_env[key]}")
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, env=run_env, cwd=cwd
        )
        logger.info(output.decode())
    except subprocess.CalledProcessError as exc:
        logger.error(f"Command execution failed: {exc.output.decode()}")
        raise SystemExit(1) from exc
