from logging import Logger
import concurrent.futures
from command_runner import run_logged_command
from thread_allocation import (
    MIN_SIPROS_OR_PERCOLATOR_THREADS,
    allocate_threads,
    available_cpu_count,
    configure_process_worker,
    effective_thread_count,
    thread_env_updates,
)

# Set native-library limits before pandas imports NumPy/OpenBLAS.  Process-pool
# workers inherit this initialized one-thread runtime and external tools receive
# their larger per-job quota in run_command().
configure_process_worker(1)

import pandas as pd


class filter:
    def __init__(self, baseNames: list[str], outputPath: str, percolatorPath: str,
                threadNumber: int, logger: Logger, ignorePCT=False, dryrun=False) -> None:
        self.percolatorPath = percolatorPath
        self.baseNames = baseNames
        self.outPutPath = outputPath
        self.logger = logger
        self.core_count = available_cpu_count()
        self.threadNumber = effective_thread_count(threadNumber, self.core_count)
        self.ignorePCT = ignorePCT
        self.dryrun = dryrun
        
    def intergrate_filtered_psms_with_feature(self):
        pass

    def run_command(self, cmd: str, threads: int):
        run_logged_command(
            cmd,
            self.logger,
            env_updates=thread_env_updates(threads),
        )
            
    def ignore_pct_in_pin(self, baseName) -> None:
        pct_columns = ['MS1IsotopicAbundances', 'MS2IsotopicAbundances', 'isotopicAbundanceDiffs']
        self.logger.info(f'Ignoring {", ".join(pct_columns)} in {baseName}_NoPCT.pin')
        pin = pd.read_csv(f'{self.outPutPath}/{baseName}/{baseName}.pin', sep='\t')
        pin.drop(columns=pct_columns, inplace=True, errors='ignore')
        pin.to_csv(f'{self.outPutPath}/{baseName}/{baseName}_NoPCT.pin', sep='\t', index=False)

    def run(self) -> None:
        # Call the feature extraction tool
        self.logger.info(f'Running Percolator: {self.percolatorPath}')
        allocation = allocate_threads(
            self.threadNumber,
            len(self.baseNames),
            minimum_threads_per_task=MIN_SIPROS_OR_PERCOLATOR_THREADS,
        )
        if allocation.worker_count == 0:
            self.logger.info('Percolator: no jobs')
            return
        minimum = min(allocation.task_threads)
        maximum = max(allocation.task_threads)
        per_job = str(minimum) if minimum == maximum else f'{minimum}-{maximum}'
        self.logger.info(
            f'Percolator: up to {allocation.worker_count} concurrent jobs, '
            f'{per_job} threads per job, {allocation.peak_threads}/{self.threadNumber} '
            f'threads in the first wave'
        )
        if not self.dryrun:
            postfix = ''
            if self.ignorePCT:
                preprocessing_allocation = allocate_threads(
                    self.threadNumber, len(self.baseNames)
                )
                with concurrent.futures.ProcessPoolExecutor(
                        max_workers=preprocessing_allocation.worker_count,
                        initializer=configure_process_worker,
                        initargs=(1,)) as executor:
                    list(executor.map(self.ignore_pct_in_pin, self.baseNames))
                postfix = '_NoPCT'
            commands = [
                f'{self.percolatorPath} --only-psms --no-terminate --num-threads {threads} '
                f'--results-psms {self.outPutPath}/{baseName}/{baseName}_target_psms.tsv '
                f'--decoy-results-psms {self.outPutPath}/{baseName}/{baseName}_decoy_psms.tsv '
                f'{self.outPutPath}/{baseName}/{baseName}{postfix}.pin'
                for baseName, threads in zip(self.baseNames, allocation.task_threads)
            ]
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=allocation.worker_count) as executor:
                list(executor.map(
                    self.run_command, commands, allocation.task_threads
                ))
