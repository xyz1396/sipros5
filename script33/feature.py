from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import re
from pathlib import Path

class feature:
    def __init__(self, baseNames: list[str], outputPath: str, scansPerFT: str, aerithFeatureExtractorPath: str,
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
            
    def get_smallest_pct_cfg(self) -> str:
        """
        Get the path of the smallest percentage .cfg file larger than 0%.
        The file name format is like SIP_C13_001.020Pct.cfg.
        The files are located in the 'configs' folder of the output path.
        """
        configs_dir = Path(self.outPutPath) / "configs"
        cfg_files = list(configs_dir.glob("*.cfg"))
        smallest_pct = None
        smallest_file = None
        # for regular search
        if len(cfg_files) == 1:
            smallest_file = str(cfg_files[0])
            return smallest_file
        # for SIP search
        for cfg_file in cfg_files:
            match = re.search(r'(\d{3}\.\d{3})Pct\.cfg$', cfg_file.name)
            if match:
                pct = float(match.group(1))  # Convert "001.020" to 1.02 for example
                if pct > 0 and (smallest_pct is None or pct < smallest_pct):
                    smallest_pct = pct
                    smallest_file = cfg_file
        if smallest_file is None:
            self.logger.error(f"No valid .cfg file found with SIP percentage > 0 in the {configs_dir}")
            exit(1)

        return str(smallest_file)

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
        
        # Get the smallest percentage .cfg file
        smallest_cfg_path = self.get_smallest_pct_cfg()
        self.logger.info(f'aerithFeatureExtractor is using config file: {smallest_cfg_path}')
        
        self.logger.info(f'Running Aerith Feature Extractor with {raw_file_parallel} processes')
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            commands = [f'{self.aerithFeatureExtractorPath} -t {self.outPutPath}/{baseName}/target \
                        -d {self.outPutPath}/{baseName}/decoy \
                        -n 8 -f {self.outPutPath}/{baseName}/ft \
                        -j {raw_file_parallel} -c {smallest_cfg_path} \
                        -p 0 -o {self.outPutPath}/{baseName}/{baseName}.pin'
                        for baseName in self.baseNames]
            if not self.dryrun:
                list(executor.map(self.run_command, commands))
