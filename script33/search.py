from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import shutil


class search:
    def __init__(self, toleranceMS1: float, toleranceMS2: float,
                 sipRange: str, step: str,
                 configGeneratorPath: str,
                 configTemplatePath: str, raxportPath: str,
                 siprosPath: str, scansPerFT2: str, fastaPath: str,
                 inputPath: str, outputPath: str, negative_control: str,
                 threadNumber: int, logger: Logger, element = "R",
                 nPrecursor = 6, dryrun=False) -> None:
        self.core_count: int = multiprocessing.cpu_count()
        self.element = element
        self.toleranceMS1 = toleranceMS1
        self.toleranceMS2 = toleranceMS2
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
        self.negative_control = negative_control
        self.threadNumber = threadNumber
        self.OMP_NUM_THREADS = 10
        self.logger = logger
        self.raw_files: list[str] = []
        self.mzml_files: list[str] = []
        self.base_names: list[str] = []
        self.base_names_of_raw: list[str] = []
        self.base_names_of_mzml: list[str] = []
        self.nPrecursor = nPrecursor
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
            
    def run_command_raxport(self, cmd):
        # Check if FT2 file exists before running raxport
        parts = cmd.split()
        try:
            f_idx = parts.index('-f')
            raw_file = parts[f_idx + 1]
            o_idx = parts.index('-o')
            ft2_dir = parts[o_idx + 1]
            if os.path.exists(ft2_dir) and any(f.endswith('.FT2') for f in os.listdir(ft2_dir)):
                self.logger.info(f"FT2 files already exist in {ft2_dir}, skipping conversion for {raw_file}")
                return
        except Exception as e:
            self.logger.error(f"Exception occurred while checking FT2 files: {e}")
            pass
        self.run_command(cmd)

    def run_command_sipros(self, cmd):
        # Try to extract output file name from command
        output_file = None
        # Look for "-o <dir>" and try to reconstruct output file name
        parts = cmd.split()
        try:
            o_idx = parts.index('-o')
            out_dir = parts[o_idx + 1]
            # Find base name and search name from command
            f_idx = parts.index('-f')
            input_file = parts[f_idx + 1]
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            c_idx = parts.index('-c')
            config_file = parts[c_idx + 1]
            if (self.element == "R"):
                search_name:str = "SE"
            else:
                search_name = os.path.splitext(os.path.basename(config_file))[0]
                search_name = search_name.replace('.', '_')
            output_file = f"{out_dir}/{base_name}.{search_name}.Spe2Pep.txt"
        except Exception as e:
            self.logger.error(f"Exception occurred while checking .Spe2Pep.txt files: {e}")
            output_file = None
            pass
        if output_file and os.path.exists(output_file) and os.path.getsize(output_file) > 500 * 1024:
            self.logger.info(f"the output file {output_file} existed, skip this search")
            return
        self.run_command(cmd)

    def reverse_fasta_sequences(self):
        self.logger.info(f'Reversing fasta sequences to {self.decoyPath}')
        if not os.path.exists(self.fastaPath):
            self.logger.error(f'Fasta file {self.fastaPath} does not exist')
            exit(1)
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
        cfgTempName = "SIP.cfg"
        toleranceMS1 = 0.01
        toleranceMS2 = 0.01
        sip_lower_bound = 0
        sip_higher_bound = 100
        sip_step = 1
        if self.toleranceMS1 != None:
            toleranceMS1 = self.toleranceMS1
        if self.toleranceMS2 != None:
            toleranceMS2 = self.toleranceMS2
        if self.sipRange != None:
            sip_lower_bound = int(self.sipRange.split('-')[0])
            sip_higher_bound = int(self.sipRange.split('-')[1])
        if self.step != None:
            sip_step = int(self.step)
        if self.element == "R":
            self.logger.info("It is regular search")
            cfgTempName = "Regular.cfg"
        else:
            self.logger.info("It is SIP search")
        configGenerator_cmd = f'{self.configGeneratorPath} \
        -i {self.configTemplatePath}/{cfgTempName} \
        -o {self.outPutPath}/configs/ \
        -t1 {toleranceMS1} -t2 {toleranceMS2} \
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
                            -j {min(10, self.core_count)} -n {self.nPrecursor}'
                self.run_command_raxport(raxport_cmd)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=min(10, self.core_count)) as executor:
                commands = [f'{self.raxportPath} -f {self.raw_files[i]} \
                            -o {self.outPutPath}/{self.base_names_of_raw[i]}/ft -n {self.nPrecursor}'
                            for i in range(len(self.raw_files))]
                list(executor.map(self.run_command_raxport, commands))

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
                                -o {self.outPutPath}/{self.base_names_of_mzml[i]}/decoy')
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            list(executor.map(self.run_command_sipros, commands))

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
        # Verify negative control files are in base_names
        if (self.negative_control != None) and (self.negative_control != ''):
            negative_control_files = self.negative_control.split(',')
            for nc_file in negative_control_files:
                nc_file = nc_file.strip()
                if nc_file.strip() not in self.base_names:
                    self.logger.error(f'Negative control file {nc_file} not found in input files')
                    self.logger.error('Please check your input files and negative control files!')
                    exit(1)
            self.logger.info(f'negative control files: {self.negative_control} verified in input files')

        for base_name in self.base_names:
            os.makedirs(f'{self.outPutPath}/{base_name}', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/ft', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/target', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/decoy', exist_ok=True)
        if not self.dryrun:
            self.convert_raw_to_ft2()
            self.sipros_search(raw_file_parallel)
