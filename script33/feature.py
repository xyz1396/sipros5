from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures


class feature:
    def __init__(self, baseNames: list[str], outputPath: str, scansPerFT: int, aerithFeatureExtractorPath: str,
                 configTemplatePath: str, threadNumber: int, logger: Logger, dryrun=False) -> None:
        self.aerithFeatureExtractorPath = aerithFeatureExtractorPath
        self.configTemplatePath = configTemplatePath
        self.baseNames = baseNames
        self.outPutPath = outputPath
        self.logger = logger
        self.OMP_NUM_THREADS = 10
        self.scansPerFT = scansPerFT
        self.threadNumber = threadNumber
        self.core_count: int = multiprocessing.cpu_count()
        self.dryrun = dryrun

    def run_command(self, cmd):
        self.logger.info(f"Running command: {cmd}")
        try:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(self.OMP_NUM_THREADS)
            # to void segmentation fault of stack overflow
            env["OMP_STACKSIZE"] = "16M"
            output = subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT, env=env)
            self.logger.info(output.decode())
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command execution failed: {e.output.decode()}")
            exit(1)

    def run(self) -> None:
        # Call the feature extraction tool
        self.logger.info(f'Runing Aerith Feature Extractor: {self.aerithFeatureExtractorPath}')
        threadNumber = self.core_count
        if self.threadNumber != None:
            threadNumber = self.threadNumber
        raw_file_parallel = threadNumber
        if self.scansPerFT != None:
            raw_file_parallel = int(threadNumber // self.OMP_NUM_THREADS)
        raw_file_parallel = min(raw_file_parallel, 10)
        self.logger.info(f'Running Aerith Feature Extractor with {raw_file_parallel} processes')
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            commands = [f'{self.aerithFeatureExtractorPath} -t {self.outPutPath}/{baseName}/target \
                        -d {self.outPutPath}/{baseName}/decoy \
                        -n 8 -f {self.outPutPath}/{baseName}/ft \
                        -j {raw_file_parallel} -c {self.configTemplatePath} \
                        -p 0 -o {self.outPutPath}/{baseName}/{baseName}.pin'
                        for baseName in self.baseNames]
            if not self.dryrun:
                list(executor.map(self.run_command, commands))
