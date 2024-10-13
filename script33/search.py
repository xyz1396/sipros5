from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import shutil

class search:
    def __init__(self, element: str, sipRange:str, step: str, 
                 configGeneratorPath: str,
                 configTemplatePath: str, raxportPath: str,  
                 siprosPath : str, scansPerFT2: str, fastaPath: str, 
                 inputPath : str, outputPath: str,
                 threadNumber: int, logger: Logger) -> None:
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
        self.DecoyPath = f'{outputPath}/Decoy.faa'
        self.inputPath = inputPath
        self.outPutPath = outputPath
        self.threadNumber = threadNumber
        self.OMP_NUM_THREADS = 10
        self.logger = logger
        self.raw_files:list[str] = []
        self.base_names:list[str] = []

    def run_command(self, cmd):
        self.logger.info(f"Running command: {cmd}")
        try:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(self.OMP_NUM_THREADS)
            output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, env=env)
            self.logger.info(output.decode())
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command execution failed: {e.output.decode()}")
            exit(1)
    
    def reverse_fasta_sequences(self):
        self.logger.info(f'Reversing fasta sequences to {self.DecoyPath}')
        with open(self.fastaPath, 'r') as fasta, \
            open(self.DecoyPath, 'w') as output:
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
            shutil.copy(f'{self.configTemplatePath}/Regular.cfg', f'{self.outPutPath}/configs/Regular.cfg')
        # if element is provided, generate config files for SIP search
        else:
            self.logger.info("It is SIP search")
            configGenerator_cmd = f'{self.configGeneratorPath} \
            -i {self.configTemplatePath}/SIP.cfg \
            -o {self.outPutPath}/configs/ \
            -e {self.element} -l {sip_lower_bound} -u {sip_higher_bound} -s {sip_step}'
            self.run_command(configGenerator_cmd)

    def convert_raw_to_ft2(self):
        self.logger.info(f'Converting raw files to FT2 files')
        if self.scansPerFT2 != None:
            scansPerFT2 = int(self.scansPerFT2)
            for base_name in self.base_names:
                raxport_cmd = f'{self.raxportPath} -f {self.inputPath}/{base_name}.raw \
                            -o {self.outPutPath}/{base_name}/ft -s {scansPerFT2} -t {min(10, self.core_count)}'
                self.run_command(raxport_cmd)            
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(10, self.core_count)) as executor:
                commands = [f'{self.raxportPath} -f {self.inputPath}/{base_name}.raw \
                            -o {self.outPutPath}/{base_name}/ft' for base_name in self.base_names]
                executor.map(self.run_command, commands)

    def sipros_search(self, raw_file_parallel: int):
        config_files = []
        for file in os.listdir(self.configsPath):
            if file.endswith(".cfg"):
                config_files.append(os.path.join(self.configsPath, file))
        # Search target database
        commands = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=raw_file_parallel) as executor:
            for base_name in self.base_names:
                for config in config_files:
                    for ft2_file in os.listdir(f'{self.outPutPath}/{base_name}/ft'):
                        if ft2_file.endswith(".FT2"):
                            commands.append(f'{self.siprosPath} -c {config} \
                                            -fasta {self.fastaPath} \
                                            -f {self.outPutPath}/{base_name}/ft/{ft2_file} \
                                            -o {self.outPutPath}/{base_name}/target')
            executor.map(self.run_command, commands)
        # Search decoy database
        commands = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=raw_file_parallel) as executor:
            for base_name in self.base_names:
                for config in config_files:
                    for ft2_file in os.listdir(f'{self.outPutPath}/{base_name}/ft'):
                        if ft2_file.endswith(".FT2"):
                            commands.append(f'{self.siprosPath} -c {config} \
                                            -fasta {self.DecoyPath} \
                                            -f {self.outPutPath}/{base_name}/ft/{ft2_file} \
                                            -o {self.outPutPath}/{base_name}/decoy')
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

        for file in os.listdir(self.inputPath):
            if file.endswith(".raw"):
                self.raw_files.append(file)
                self.base_names.append(os.path.splitext(file)[0])
        if len(self.raw_files) == 0:
            self.logger.error(f'No raw files found in {self.inputPath}')
            exit(1)
        self.logger.info(f'Raw files: {self.raw_files}')    
        
        for base_name in self.base_names:
            os.makedirs(f'{self.outPutPath}/{base_name}', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/ft', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/target', exist_ok=True)
            os.makedirs(f'{self.outPutPath}/{base_name}/decoy', exist_ok=True)

        self.convert_raw_to_ft2() 

        self.sipros_search(raw_file_parallel)           

        

