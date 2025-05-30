from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import pandas as pd

class filter:
    def __init__(self, baseNames: list[str], outputPath: str, percolatorPath: str,
                threadNumber: int, logger: Logger, ignorePCT=False, dryrun=False) -> None:
        self.percolatorPath = percolatorPath
        self.baseNames = baseNames
        self.outPutPath = outputPath
        self.logger = logger
        self.OMP_NUM_THREADS = 10
        self.threadNumber = threadNumber
        self.core_count: int = multiprocessing.cpu_count()
        self.ignorePCT = ignorePCT
        self.dryrun = dryrun
        
    def intergrate_filtered_psms_with_feature(self):
        pass

    def run_command(self, cmd):
        self.logger.info(f"Running command: {cmd}")
        try:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(self.OMP_NUM_THREADS)
            output = subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT, env=env)
            self.logger.info(output.decode())
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command execution failed: {e.output.decode()}")
            exit(1)
            
    def ignore_pct_in_pin(self, baseName) -> None:
        self.logger.info(f'Ignoring MS1IsotopicAbundances MS2IsotopicAbundances in {baseName}_NoPCT.pin')
        pin = pd.read_csv(f'{self.outPutPath}/{baseName}/{baseName}.pin', sep='\t')
        pin.drop(columns=['MS1IsotopicAbundances', 'MS2IsotopicAbundances'], inplace=True)
        pin.to_csv(f'{self.outPutPath}/{baseName}/{baseName}_NoPCT.pin', sep='\t', index=False)

    def run(self) -> None:
        # Call the feature extraction tool
        self.logger.info(f'Runing percolator: {self.percolatorPath}')
        threadNumber = self.core_count
        if self.threadNumber != None:
            threadNumber = self.threadNumber
        raw_file_parallel = int(threadNumber // self.OMP_NUM_THREADS)
        self.logger.info(f'Running percolator with {raw_file_parallel} processes')
        # if not self.dryrun:
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            postfix = ''
            if self.ignorePCT:
                list(executor.map(self.ignore_pct_in_pin, self.baseNames))
                postfix = '_NoPCT'
            commands = [f'{self.percolatorPath} --only-psms --no-terminate --num-threads {self.OMP_NUM_THREADS} \
                    --results-psms {self.outPutPath}/{baseName}/{baseName}_target_psms.tsv \
                    --decoy-results-psms {self.outPutPath}/{baseName}/{baseName}_decoy_psms.tsv \
                    {self.outPutPath}/{baseName}/{baseName}{postfix}.pin' for baseName in self.baseNames]
            list(executor.map(self.run_command, commands))
