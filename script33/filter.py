from logging import Logger
import shlex

from command_runner import run_logged_command
from thread_allocation import (
    available_cpu_count,
    configure_process_worker,
    effective_thread_count,
    thread_env_updates,
)

# Keep NumPy/OpenBLAS single-threaded if another workflow module imports
# pandas after this one. Aerith itself receives the complete native thread team.
configure_process_worker(1)


class filter:
    """Run one cross-sample Aerith target/decoy filtering job."""

    def __init__(self, baseNames: list[str], outputPath: str, aerithPath: str,
                 threadNumber: int, logger: Logger, decoyPrefix: str = "Decoy_",
                 ignorePCT: bool = False, dryrun: bool = False,
                 spectraPaths: list[str] | None = None) -> None:
        self.aerithPath = aerithPath
        self.baseNames = baseNames
        self.outPutPath = outputPath
        self.logger = logger
        self.core_count = available_cpu_count()
        self.threadNumber = effective_thread_count(threadNumber, self.core_count)
        self.decoyPrefix = decoyPrefix
        self.ignorePCT = ignorePCT
        self.dryrun = dryrun
        self.spectraPaths = [] if spectraPaths is None else list(spectraPaths)
        if self.spectraPaths and len(self.spectraPaths) != len(self.baseNames):
            raise ValueError("Aerith requires one spectra path per sample")

    def command(self) -> str:
        arguments = [
            self.aerithPath,
            "--database", f"{self.outPutPath}/targetDecoy.faa",
            "--decoy-prefix", self.decoyPrefix,
        ]
        if self.ignorePCT:
            arguments.append("--ignore-pct")
        for index, baseName in enumerate(self.baseNames):
            sample = f"{self.outPutPath}/{baseName}/{baseName}"
            arguments.extend([
                "--target-pin", f"{sample}_target.pin",
                "--decoy-pin", f"{sample}_decoy.pin",
                "--output-prefix", sample,
            ])
            if self.spectraPaths:
                arguments.extend(["--spectra", self.spectraPaths[index]])
        return shlex.join(arguments)

    def run(self) -> None:
        if not self.baseNames:
            self.logger.info("Aerith: no jobs")
            return
        command = self.command()
        self.logger.info(
            f"Running Aerith cross-sample SVM+RT filtering with "
            f"{self.threadNumber} CPU cores"
        )
        self.logger.info(command)
        if self.dryrun:
            return
        environment = thread_env_updates(self.threadNumber)
        # LibTorch gives MKL_NUM_THREADS precedence over OMP_NUM_THREADS when
        # initializing its OpenMP pool. Aerith runs as one process, so give
        # spectrum prediction the full CPU team allocated to this job.
        environment["MKL_NUM_THREADS"] = str(self.threadNumber)
        run_logged_command(
            command,
            self.logger,
            env_updates=environment,
            cpu_cores=self.threadNumber,
        )
