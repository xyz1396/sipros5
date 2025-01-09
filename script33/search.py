from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import shutil


class search:
    def __init__(self, element: str, sipRange: str, step: str,
                 configGeneratorPath: str,
                 configTemplatePath: str, raxportPath: str,
                 siprosPath: str, scansPerFT2: str, fastaPath: str,
                 inputPath: str, outputPath: str,
                 threadNumber: int, logger: Logger, dryrun=False) -> None:
        self.core_count: int = multiprocessing.cpu_count()
        self.element = element
        self.sipRange = sipRange
        self.step = step
        self.configTemplatePath = configTemplatePath
        self.configsPath = f'{outputPath}/configs'
        self.configGeneratorPath = configGeneratorPath
        self.raxportPath = raxportPath
        self.scansPerFT2 = scansPerFT2
        self.siprosPath = siprosPath
        self.fastaPath = fastaPath
        self.decoyPath = f'{outputPath}/decoy.faa'
        self.inputPath = inputPath
        self.outPutPath = outputPath
        self.threadNumber = threadNumber
        self.OMP_NUM_THREADS = 10
        self.logger = logger
        self.raw_files: list[str] = []
        self.mzml_files: list[str] = []
        self.base_names: list[str] = []
        self.base_names_of_raw: list[str] = []
        self.base_names_of_mzml: list[str] = []
        self.dryrun = dryrun

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

    def reverse_fasta_sequences(self):
        self.logger.info(f'Reversing fasta sequences to {self.decoyPath}')
        with open(self.fastaPath, 'r') as fasta, \
                open(self.decoyPath, 'w') as output:
            sequence = ''
            header = ''
            for line in fasta:
                if line.startswith('>'):
                    if header:
                        output.write(header + '\n' + sequence[::-1] + '\n')
                    header = '>Decoy_' + line[1:].strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            if header:  # write the last sequence
                output.write(header + '\n' + sequence[::-1] + '\n')

    def generateConfigs(self):
        self.logger.info(f'Generating config files to {self.configsPath}')
        if not os.path.exists(self.configsPath):
            os.makedirs(self.configsPath)
        sip_lower_bound = 0
        sip_higher_bound = 100
        sip_step = 1
        if self.sipRange != None:
            sip_lower_bound = int(self.sipRange.split('-')[0])
            sip_higher_bound = int(self.sipRange.split('-')[1])
        if self.step != None:
            sip_step = int(self.step)
        if self.element == None:
            self.logger.info("It is regular search")
            shutil.copy(f'{self.configTemplatePath}/Regular.cfg',
                        f'{self.outPutPath}/configs/Regular.cfg')
        # if element is provided, generate config files for SIP search
        else:
            self.logger.info("It is SIP search")
            configGenerator_cmd = f'{self.configGeneratorPath} \
            -i {self.configTemplatePath}/SIP.cfg \
            -o {self.outPutPath}/configs/ \
            -e {self.element} -l {sip_lower_bound} -u {sip_higher_bound} -s {sip_step}'
            self.run_command(configGenerator_cmd)

    def getInputFiles(self):
        files = []
        if os.path.isdir(self.inputPath):
            self.logger.info(f'{self.inputPath} is a directory')
            for file in os.listdir(self.inputPath):
                files.append(f'{self.inputPath}/{file}')
        else:
            self.logger.info(f'{self.inputPath} is a file list')
            files = self.inputPath.split(',')
        for file in files:
            if not os.path.exists(file):
                self.logger.error(f'{file} does not exist')
                exit(1)
            if file.endswith(".raw"):
                self.raw_files.append(file)
                # let base_names_of_raw match the raw files path
                raw_base = os.path.splitext(os.path.basename(file))[0]
                self.base_names_of_raw.append(raw_base)
                self.base_names.append(raw_base)
            if file.endswith(".mzml"):
                self.mzml_files.append(file)
                # let base_names_of_raw match the mzml files path
                mzml_base = os.path.splitext(os.path.basename(file))[0]
                self.base_names.append(mzml_base)
                self.base_names_of_mzml.append(mzml_base)
        if len(self.raw_files) == 0 and len(self.mzml_files) == 0:
            self.logger.error(
                f'No raw or mzml files found in {self.inputPath}')
            exit(1)
        self.logger.info(f'raw files: {self.raw_files}')
        self.logger.info(f'mzml files: {self.mzml_files}')

    def convert_raw_to_ft2(self):
        self.logger.info(f'Converting raw files to FT2 files')
        # split FT2 files
        if self.scansPerFT2 != None:
            scansPerFT2 = int(self.scansPerFT2)
            for i in range(len(self.raw_files)):
                raxport_cmd = f'{self.raxportPath} -f {self.raw_files[i]} \
                            -o {self.outPutPath}/{self.base_names_of_raw[i]}/ft -s {scansPerFT2} \
                            -j {min(10, self.core_count)}'
                self.run_command(raxport_cmd)
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(10, self.core_count)) as executor:
                commands = [f'{self.raxportPath} -f {self.raw_files[i]} \
                            -o {self.outPutPath}/{self.base_names_of_raw[i]}/ft'
                            for i in range(len(self.raw_files))]
                executor.map(self.run_command, commands)

    def sipros_search(self, raw_file_parallel: int):
        config_files = []
        for file in os.listdir(self.configsPath):
            if file.endswith(".cfg"):
                config_files.append(os.path.join(self.configsPath, file))
        commands = []
        # Search target database
        for config in config_files:
            for base_name in self.base_names_of_raw:
                for ft2_file in os.listdir(f'{self.outPutPath}/{base_name}/ft'):
                    if ft2_file.endswith(".FT2"):
                        commands.append(f'{self.siprosPath} -c {config} \
                                        -fasta {self.fastaPath} \
                                        -f {self.outPutPath}/{base_name}/ft/{ft2_file} \
                                        -o {self.outPutPath}/{base_name}/target')
            for i in range(len(self.mzml_files)):
                commands.append(f'{self.siprosPath} -c {config} \
                                -fasta {self.fastaPath} \
                                -f {self.mzml_files[i]} \
                                -o {self.outPutPath}/{self.base_names_of_mzml[i]}/target')
        # Search decoy database
        for config in config_files:
            for base_name in self.base_names_of_raw:
                for ft2_file in os.listdir(f'{self.outPutPath}/{base_name}/ft'):
                    if ft2_file.endswith(".FT2"):
                        commands.append(f'{self.siprosPath} -c {config} \
                                        -fasta {self.decoyPath} \
                                        -f {self.outPutPath}/{base_name}/ft/{ft2_file} \
                                        -o {self.outPutPath}/{base_name}/decoy')
            for i in range(len(self.mzml_files)):
                commands.append(f'{self.siprosPath} -c {config} \
                                -fasta {self.decoyPath} \
                                -f {self.mzml_files[i]} \
                                -o {self.outPutPath}/{self.base_names_of_mzml[i]}/target')
        with concurrent.futures.ThreadPoolExecutor(max_workers=raw_file_parallel) as executor:
            executor.map(self.run_command, commands)

    def run(self) -> None:
        self.reverse_fasta_sequences()
        self.generateConfigs()
        self.logger.info(f'Number of CPU cores: {self.core_count}')
        threadNumber = self.core_count
        if self.threadNumber != None:
            threadNumber = self.threadNumber
        self.logger.info(f'Setted max thread numbers: {threadNumber}')
        raw_file_parallel = int(threadNumber // self.OMP_NUM_THREADS)
        self.getInputFiles()

        for base_name in self.base_names:
            os.makedirs(f'{self.outPutPath}/{base_name}', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/ft', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/target', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/decoy', exist_ok=True)
        if not self.dryrun:
            self.convert_raw_to_ft2()
            self.sipros_search(raw_file_parallel)
