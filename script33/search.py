from logging import Logger
import os
import concurrent.futures
import shutil
import shlex
from pathlib import Path
from command_runner import run_logged_command
from thread_allocation import (
    MIN_SIPROS_OR_PERCOLATOR_THREADS,
    ThreadAllocation,
    allocate_threads,
    available_cpu_count,
    effective_thread_count,
    thread_env_updates,
)


class search:
    def __init__(self, toleranceMS1: float, toleranceMS2: float,
                 sipRange: str, step: str,
                 configTemplatePath: str, raxportPath: str,
                 siprosPath: str, fastaPath: str,
                 inputPath: str, outputPath: str, negative_control: str,
                 threadNumber: int, logger: Logger, element="R",
                 nPrecursor=6, dryrun=False, psmTsv: str | None = None,
                 unlabeledInput: str | None = None, spectraDir: str | None = None,
                 topPsmsPerScan: int = 8) -> None:
        self.core_count = available_cpu_count()
        self.element = element
        self.toleranceMS1 = toleranceMS1
        self.toleranceMS2 = toleranceMS2
        self.sipRange = sipRange
        self.step = step
        self.configTemplatePath = configTemplatePath
        self.raxportPath = raxportPath
        self.siprosPath = siprosPath
        self.fastaPath = fastaPath
        self.decoyPath = f'{outputPath}/decoy.faa'
        self.inputPath = inputPath
        self.outPutPath = outputPath
        self.negative_control = negative_control
        self.threadNumber = effective_thread_count(threadNumber, self.core_count)
        self.logger = logger
        self.raw_files: list[str] = []
        self.hdf5_input_files: list[str] = []
        self.base_names: list[str] = []
        self.base_names_of_raw: list[str] = []
        self.base_names_of_hdf5: list[str] = []
        self.hdf5_paths: dict[str, str] = {}
        self.nPrecursor = nPrecursor
        self.dryrun = dryrun
        self.psmTsv = psmTsv
        self.unlabeledInput = unlabeledInput
        self.spectraDir = spectraDir
        self.topPsmsPerScan = topPsmsPerScan
        self.generatedSpectraDir = f'{outputPath}/spectra'
        self.decoyPrefix = 'DECOY_' if (self.psmTsv or self.unlabeledInput or self.spectraDir) else 'Decoy_'

    def q(self, value: str | Path) -> str:
        return shlex.quote(str(value))

    def run_command(self, cmd: str, env: dict[str, str] | None = None,
                    threads: int | None = None):
        thread_count = self.threadNumber if threads is None else threads
        run_logged_command(
            cmd,
            self.logger,
            env=env,
            env_updates=thread_env_updates(thread_count),
            cpu_cores=thread_count,
        )

    def log_thread_allocation(self, phase: str, allocation: ThreadAllocation) -> None:
        if allocation.worker_count == 0:
            self.logger.info(f'{phase}: no jobs')
            return
        minimum = min(allocation.task_threads)
        maximum = max(allocation.task_threads)
        per_job = str(minimum) if minimum == maximum else f'{minimum}-{maximum}'
        unit = 'core' if minimum == maximum == 1 else 'cores'
        self.logger.info(
            f'{phase}: up to {allocation.worker_count} concurrent processes; '
            f'{per_job} CPU {unit} per process; '
            f'{allocation.peak_threads}/{self.threadNumber} cores allocated at peak'
        )

    def sample_base_name(self, path: str) -> str:
        name = Path(path).name
        lower = name.lower()
        if lower.endswith('.d.zip'):
            return name[:-6]
        return Path(name).stem

    def is_raw_input(self, path: str) -> bool:
        lower = path.lower()
        return lower.endswith('.raw') or lower.endswith('.d') or lower.endswith('.d.zip')

    def is_hdf5_input(self, path: str) -> bool:
        lower = path.lower()
        return lower.endswith('.h5') or lower.endswith('.hdf5')

    def expected_hdf5_path(self, base_name: str) -> str:
        return f'{self.outPutPath}/{base_name}/{base_name}.h5'

    def complete_hdf5_exists(self, path: str) -> bool:
        return os.path.exists(path) and os.path.getsize(path) >= 1024 * 1024

    def run_command_raxport(self, raw_file: str, hdf5_dir: str,
                            expected_hdf5: str, threads: int):
        if os.path.exists(expected_hdf5):
            hdf5_size = os.path.getsize(expected_hdf5)
            if self.complete_hdf5_exists(expected_hdf5):
                self.logger.info(f'HDF5 file already exists, skipping conversion: {expected_hdf5}')
                return
            self.logger.warning(f'Removing incomplete HDF5 conversion output: {expected_hdf5} ({hdf5_size} bytes)')
            os.remove(expected_hdf5)
        os.makedirs(hdf5_dir, exist_ok=True)
        env = os.environ.copy()
        raxport_heap_limit = (
            os.environ.get("SIPROS_RAXPORT_GC_HEAP_LIMIT")
            or os.environ.get("DOTNET_GCHeapHardLimit")
            or str(128 * 1024 * 1024)
        )
        env["DOTNET_GCHeapHardLimit"] = raxport_heap_limit
        self.logger.debug(f"Raxport heap limit: {raxport_heap_limit} bytes")
        cmd = (f'{self.q(self.raxportPath)} -f {self.q(raw_file)} -o {self.q(hdf5_dir)} '
               f'--format hdf5 -n {self.nPrecursor}')
        self.run_command(cmd, env, threads)
        if not os.path.exists(expected_hdf5):
            raise FileNotFoundError(f'Raxport did not create expected HDF5 file: {expected_hdf5}')

    def run_command_sipros(self, cmd: str, threads: int | None = None):
        self.run_command(cmd, threads=threads)

    def reverse_fasta_sequences(self):
        self.logger.info(f'Reversing fasta sequences to {self.decoyPath}')
        if not os.path.exists(self.fastaPath):
            self.logger.error(f'Fasta file {self.fastaPath} does not exist')
            raise SystemExit(1)
        with open(self.fastaPath, 'r') as fasta, open(self.decoyPath, 'w') as output:
            sequence = ''
            header = ''
            for line in fasta:
                if line.startswith('>'):
                    if header:
                        output.write(header + '\n' + sequence[::-1] + '\n')
                    header = '>' + self.decoyPrefix + line[1:].strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            if header:
                output.write(header + '\n' + sequence[::-1] + '\n')

    def sip_element_symbol(self) -> str:
        if self.element == "R" or not self.element:
            return "C"
        return self.element[0].upper()

    def sip_isotope_number(self) -> str:
        if self.element == "R" or not self.element:
            return "13"
        supported = {
            "C13": "13",
            "H2": "2",
            "N15": "15",
            "O18": "18",
            "S34": "34",
        }
        try:
            return supported[self.element]
        except KeyError as exc:
            raise ValueError(
                f"Unsupported SIP isotope {self.element}; "
                "use C13, H2, N15, O18, or S34") from exc

    def update_config_line(self, line: str, replacements: dict[str, str]) -> str:
        stripped = line.lstrip()
        active = stripped.split("#", 1)[0]
        if "=" not in active:
            return line
        line_key = active.split("=", 1)[0].strip()
        if line_key not in replacements:
            return line
        indent = line[:len(line) - len(stripped)]
        suffix = ""
        if "#" in stripped:
            suffix = "  #" + stripped.split("#", 1)[1].rstrip("\n")
        return f"{indent}{line_key} = {replacements[line_key]}{suffix}\n"

    def write_workflow_config(self) -> str:
        cfgTempName = "Regular.cfg" if self.element == "R" else "SIP.cfg"
        template_path = Path(self.configTemplatePath) / cfgTempName
        if not template_path.exists():
            self.logger.error(f"Config template does not exist: {template_path}")
            raise SystemExit(1)
        output_path = Path(self.outPutPath) / cfgTempName
        self.logger.info(f"Writing workflow config to {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        toleranceMS1 = self.toleranceMS1 if self.toleranceMS1 is not None else 0.01
        toleranceMS2 = self.toleranceMS2 if self.toleranceMS2 is not None else 0.01
        search_type = "Regular" if self.element == "R" else "SIP"
        search_name = "SE" if self.element == "R" else "SIP"
        replacements = {
            "Search_Name": search_name,
            "Search_Type": search_type,
            "FASTA_Database": self.fastaPath,
            "Search_Mass_Tolerance_Parent_Ion": str(toleranceMS1),
            "Mass_Tolerance_Fragment_Ions": str(toleranceMS2),
            "SIP_Element": self.sip_element_symbol(),
            "SIP_Element_Isotope": self.sip_isotope_number(),
        }
        with open(template_path, "r") as source, open(output_path, "w") as output:
            for line in source:
                output.write(self.update_config_line(line, replacements))
        self.logger.info(f"Wrote config file: {output_path}")
        return str(output_path)

    def input_entries(self, input_path: str) -> list[str]:
        path = Path(input_path)
        if self.is_raw_input(input_path) or self.is_hdf5_input(input_path):
            return [input_path]
        if path.is_dir():
            self.logger.info(f'{input_path} is a directory')
            return [str(path / name) for name in sorted(os.listdir(path))]
        self.logger.info(f'{input_path} is a file list')
        return [item.strip() for item in input_path.split(',') if item.strip()]

    def getInputFiles(self):
        self.raw_files.clear()
        self.hdf5_input_files.clear()
        self.base_names.clear()
        self.base_names_of_raw.clear()
        self.base_names_of_hdf5.clear()
        for file in self.input_entries(self.inputPath):
            real_file = os.path.realpath(file)
            if real_file != file:
                self.logger.info(f'{file} is not an absolute path, resolved to {real_file}')
            if not os.path.exists(real_file):
                self.logger.error(f'{real_file} does not exist')
                raise SystemExit(1)
            base = self.sample_base_name(real_file)
            if self.is_raw_input(real_file):
                self.raw_files.append(real_file)
                self.base_names_of_raw.append(base)
            elif self.is_hdf5_input(real_file):
                self.hdf5_input_files.append(real_file)
                self.base_names_of_hdf5.append(base)
            else:
                self.logger.warning(f'Skipping unsupported input: {real_file}')
                continue
            self.base_names.append(base)
        if len(self.raw_files) == 0 and len(self.hdf5_input_files) == 0:
            self.logger.error(f'No raw/.d/.d.zip or Raxport HDF5 files found in {self.inputPath}')
            raise SystemExit(1)
        self.logger.info(f'raw files: {self.raw_files}')
        self.logger.info(f'HDF5 files: {self.hdf5_input_files}')

    def create_sample_directories(self):
        for base_name in self.base_names:
            os.makedirs(f'{self.outPutPath}/{base_name}', exist_ok=True)

    def prepare_hdf5_inputs(self):
        self.logger.info('Preparing Raxport HDF5 scan inputs')
        self.hdf5_paths.clear()
        for hdf5_file, base in zip(self.hdf5_input_files, self.base_names_of_hdf5):
            expected = self.expected_hdf5_path(base)
            os.makedirs(os.path.dirname(expected), exist_ok=True)
            if os.path.realpath(hdf5_file) != os.path.realpath(expected):
                if not os.path.exists(expected):
                    self.logger.info(f'Staging HDF5 input {hdf5_file} to {expected}')
                    if not self.dryrun:
                        shutil.copy2(hdf5_file, expected)
            self.hdf5_paths[base] = expected
        if self.raw_files:
            jobs = []
            for raw_file, base in zip(self.raw_files, self.base_names_of_raw):
                expected = self.expected_hdf5_path(base)
                hdf5_dir = os.path.dirname(expected)
                if self.complete_hdf5_exists(expected):
                    self.logger.info(
                        f'HDF5 file already exists, skipping conversion: {expected}'
                    )
                    self.hdf5_paths[base] = expected
                    continue
                jobs.append((raw_file, hdf5_dir, expected, base))
            allocation = allocate_threads(self.threadNumber, len(jobs))
            self.log_thread_allocation('Raxport conversion', allocation)
            if self.dryrun:
                for _, _, expected, base in jobs:
                    self.hdf5_paths[base] = expected
            elif allocation.worker_count > 0:
                with concurrent.futures.ThreadPoolExecutor(
                        max_workers=allocation.worker_count) as executor:
                    futures = [
                        executor.submit(
                            self.run_command_raxport, raw, outdir, expected, threads
                        )
                        for (raw, outdir, expected, _), threads
                        in zip(jobs, allocation.task_threads)
                    ]
                    for future in concurrent.futures.as_completed(futures):
                        future.result()
                for _, _, expected, base in jobs:
                    self.hdf5_paths.setdefault(base, expected)
        self.logger.info(f'HDF5 scan files: {self.hdf5_paths}')

    def validate_negative_controls(self):
        if self.negative_control is not None and self.negative_control != '':
            negative_control_files = [item.strip() for item in self.negative_control.split(',') if item.strip()]
            for nc_file in negative_control_files:
                if nc_file not in self.base_names:
                    self.logger.error(f'Negative control file {nc_file} not found in input files')
                    raise SystemExit(1)
            self.logger.info(f'negative control files: {self.negative_control} verified in input files')

    def direct_sip_args(self) -> str:
        if self.element == 'R':
            return ''
        sip_range = self.sipRange if self.sipRange is not None else '0-100'
        sip_step = self.step if self.step is not None else '1'
        return (f' -a {self.q(self.element)} -b {self.q(sip_range)}'
                f' -s {self.q(sip_step)}')

    def pin_score_value(self, row: list[str], score_idx: int) -> float:
        try:
            return float(row[score_idx])
        except (ValueError, IndexError):
            return float("-inf")

    def pin_spec_id_with_rank(self, spec_id: str, rank: int) -> str:
        if "." not in spec_id:
            return spec_id
        prefix, _ = spec_id.rsplit(".", 1)
        return f"{prefix}.{rank}"

    def merge_pin_files(self, target_pin: str, decoy_pin: str, final_pin: str):
        os.makedirs(os.path.dirname(final_pin), exist_ok=True)
        header: str | None = None
        rows: list[list[str]] = []
        for pin_path in (target_pin, decoy_pin):
            if not os.path.exists(pin_path):
                raise FileNotFoundError(f"Expected direct search PIN was not created: {pin_path}")
            with open(pin_path, "r") as pin_file:
                current_header = pin_file.readline()
                if not current_header:
                    raise ValueError(f"PIN file is empty: {pin_path}")
                if header is None:
                    header = current_header
                elif current_header != header:
                    raise ValueError(f"PIN header mismatch while merging {pin_path}")
                for line in pin_file:
                    if line.strip():
                        rows.append(line.rstrip("\n").split("\t"))
        if header is None:
            raise ValueError("No PIN rows were available to merge")

        columns = header.rstrip("\n").split("\t")
        column_index = {name: idx for idx, name in enumerate(columns)}
        required_columns = ["SpecId", "ScanNr", "ranks", "WDPscores"]
        missing_columns = [name for name in required_columns if name not in column_index]
        if missing_columns:
            raise ValueError(f"Merged PIN is missing required columns: {missing_columns}")

        spec_idx = column_index["SpecId"]
        scan_idx = column_index["ScanNr"]
        rank_idx = column_index["ranks"]
        score_idx = column_index["WDPscores"]
        diff_idx = column_index.get("diffScores")
        grouped_rows: dict[str, list[list[str]]] = {}
        for row in rows:
            if len(row) != len(columns):
                raise ValueError(f"PIN row has {len(row)} fields but expected {len(columns)}: {row}")
            grouped_rows.setdefault(row[scan_idx], []).append(row)

        def scan_sort_key(scan: str):
            try:
                return (0, int(scan))
            except ValueError:
                return (1, scan)

        merged_rows: list[list[str]] = []
        for scan in sorted(grouped_rows, key=scan_sort_key):
            scan_rows = grouped_rows[scan]
            scan_rows.sort(key=lambda row: self.pin_score_value(row, score_idx), reverse=True)
            top_score = self.pin_score_value(scan_rows[0], score_idx) if scan_rows else 0.0
            for rank, row in enumerate(scan_rows, start=1):
                row[rank_idx] = str(rank)
                row[spec_idx] = self.pin_spec_id_with_rank(row[spec_idx], rank)
                if diff_idx is not None:
                    row[diff_idx] = f"{top_score - self.pin_score_value(row, score_idx):.10g}"
                merged_rows.append(row)

        with open(final_pin, "w") as merged:
            merged.write(header)
            for row in merged_rows:
                merged.write("\t".join(row) + "\n")

    def sipros_search(self, config: str):
        sip_args = self.direct_sip_args()
        direct_top_psms_per_scan = self.topPsmsPerScan
        commands: list[str] = []
        merge_jobs: list[tuple[str, str, str]] = []
        for base_name in self.base_names:
            hdf5_path = self.hdf5_paths[base_name]
            sample_dir = f'{self.outPutPath}/{base_name}'
            target_pin_name = f'{base_name}_target.pin'
            decoy_pin_name = f'{base_name}_decoy.pin'
            target_pin = f'{sample_dir}/{target_pin_name}'
            decoy_pin = f'{sample_dir}/{decoy_pin_name}'
            final_pin = f'{sample_dir}/{base_name}.pin'
            commands.append(
                f'{self.q(self.siprosPath)} search-fasta -c {self.q(config)} -fasta {self.q(self.fastaPath)} '
                f'-f {self.q(hdf5_path)} -o {self.q(sample_dir)} --pin-output {self.q(target_pin_name)}{sip_args} '
                f'--pin-label 1 --top-psms-per-scan {direct_top_psms_per_scan}'
            )
            commands.append(
                f'{self.q(self.siprosPath)} search-fasta -c {self.q(config)} -fasta {self.q(self.decoyPath)} '
                f'-f {self.q(hdf5_path)} -o {self.q(sample_dir)} --pin-output {self.q(decoy_pin_name)}{sip_args} '
                f'--pin-label -1 --top-psms-per-scan {direct_top_psms_per_scan}'
            )
            merge_jobs.append((target_pin, decoy_pin, final_pin))
        allocation = allocate_threads(
            self.threadNumber,
            len(commands),
            minimum_threads_per_task=MIN_SIPROS_OR_PERCOLATOR_THREADS,
        )
        self.log_thread_allocation('Sipros FASTA search', allocation)
        if allocation.worker_count > 0:
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=allocation.worker_count) as executor:
                futures = [
                    executor.submit(
                        self.run_command_sipros, cmd, threads
                    )
                    for cmd, threads
                    in zip(commands, allocation.task_threads)
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()
        for target_pin, decoy_pin, final_pin in merge_jobs:
            self.merge_pin_files(target_pin, decoy_pin, final_pin)

    def resolve_or_convert_unlabeled_hdf5(self) -> str:
        if self.unlabeledInput is None or self.unlabeledInput == '':
            self.logger.error('--unlabeled-input is required when --spectra-dir is not provided')
            raise SystemExit(1)
        entry = self.input_entries(self.unlabeledInput)[0]
        real_entry = os.path.realpath(entry)
        if not os.path.exists(real_entry):
            self.logger.error(f'Unlabeled input does not exist: {real_entry}')
            raise SystemExit(1)
        if self.is_hdf5_input(real_entry):
            return real_entry
        if not self.is_raw_input(real_entry):
            self.logger.error(f'Unsupported unlabeled input {real_entry}; use raw/.d/.d.zip or Raxport HDF5')
            raise SystemExit(1)
        base = self.sample_base_name(real_entry)
        out_dir = f'{self.outPutPath}/unlabeled_hdf5'
        expected = f'{out_dir}/{base}.h5'
        if not self.dryrun:
            self.run_command_raxport(
                real_entry, out_dir, expected, self.threadNumber
            )
        return expected

    def generate_or_reuse_spectra_library(self, config: str) -> str:
        if self.spectraDir:
            spectra_dir = os.path.realpath(self.spectraDir)
            if not os.path.isdir(spectra_dir):
                self.logger.error(f'--spectra-dir does not exist or is not a directory: {spectra_dir}')
                raise SystemExit(1)
            h5_files = list(Path(spectra_dir).glob('*.h5')) + list(Path(spectra_dir).glob('*.hdf5'))
            if not h5_files:
                self.logger.error(f'No HDF5 spectra library files found in --spectra-dir: {spectra_dir}')
                raise SystemExit(1)
            self.logger.info(f'Reusing HDF5 spectra library from {spectra_dir}')
            return spectra_dir
        if not self.psmTsv:
            self.logger.error('--psm-tsv is required to generate a spectra library')
            raise SystemExit(1)
        unlabeled_hdf5 = self.resolve_or_convert_unlabeled_hdf5()
        os.makedirs(self.generatedSpectraDir, exist_ok=True)
        sip_range = self.sipRange if self.sipRange is not None else '0-100'
        sip_step = self.step if self.step is not None else '1'
        cmd = (f'{self.q(self.siprosPath)} experimental-spectra -c {self.q(config)} '
               f'-i {self.q(self.psmTsv)} -f {self.q(unlabeled_hdf5)} '
               f'-o {self.q(self.generatedSpectraDir)} -a {self.element} '
               f'-b {sip_range} -s {sip_step} --decoy -t {self.threadNumber}')
        if not self.dryrun:
            self.run_command(cmd, threads=self.threadNumber)
            spectra_files = list(Path(self.generatedSpectraDir).glob('*.h5')) + list(Path(self.generatedSpectraDir).glob('*.hdf5'))
            if not spectra_files:
                self.logger.error(f'experimental-spectra did not create HDF5 spectra files in {self.generatedSpectraDir}')
                raise SystemExit(1)
            non_hdf5 = [path for path in Path(self.generatedSpectraDir).iterdir()
                        if path.is_file() and path.suffix.lower() not in {'.h5', '.hdf5'}]
            if non_hdf5:
                self.logger.error(f'Unexpected non-HDF5 spectra intermediates generated: {non_hdf5}')
                raise SystemExit(1)
        return self.generatedSpectraDir

    def search_spectra_samples(self, config: str, spectra_dir: str):
        commands: list[tuple[str, int]] = []
        samples: list[tuple[str, str, str]] = []
        for base_name in self.base_names:
            hdf5_path = self.hdf5_paths[base_name]
            sample_dir = f'{self.outPutPath}/{base_name}'
            os.makedirs(sample_dir, exist_ok=True)
            samples.append((base_name, hdf5_path, sample_dir))
        allocation = allocate_threads(
            self.threadNumber,
            len(samples),
            minimum_threads_per_task=MIN_SIPROS_OR_PERCOLATOR_THREADS,
        )
        self.log_thread_allocation('Sipros spectra search', allocation)
        if allocation.worker_count == 0:
            return
        for (_, hdf5_path, sample_dir), threads in zip(
                samples, allocation.task_threads):
            cmd = (f'{self.q(self.siprosPath)} search-spectra -f {self.q(hdf5_path)} '
                   f'-c {self.q(config)} -h5 {self.q(spectra_dir)} -o {self.q(sample_dir)} '
                   f'-t {threads} --tolerance-ms1 {self.toleranceMS1} --tolerance-ms1-unit da '
                   f'--tolerance-ms2 {self.toleranceMS2} --tolerance-ms2-unit da --top-psms-per-scan {self.topPsmsPerScan}')
            commands.append((cmd, threads))
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=allocation.worker_count) as executor:
            futures = [
                executor.submit(self.run_command_sipros, cmd, threads)
                for cmd, threads in commands
            ]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    def run_search_spectra(self) -> None:
        if self.element == 'R':
            self.logger.error('search-spectra mode requires a SIP element such as C13')
            raise SystemExit(1)
        self.reverse_fasta_sequences()
        config = self.write_workflow_config()
        self.logger.info(
            f'Workflow CPU allocation: {self.threadNumber} cores '
            f'({self.core_count} available)'
        )
        self.getInputFiles()
        self.validate_negative_controls()
        self.create_sample_directories()
        if not self.dryrun:
            self.prepare_hdf5_inputs()
            spectra_dir = self.generate_or_reuse_spectra_library(config)
            self.search_spectra_samples(config, spectra_dir)

    def run(self) -> None:
        self.reverse_fasta_sequences()
        config = self.write_workflow_config()
        self.logger.info(
            f'Workflow CPU allocation: {self.threadNumber} cores '
            f'({self.core_count} available)'
        )
        self.getInputFiles()
        self.validate_negative_controls()
        self.create_sample_directories()
        if not self.dryrun:
            self.prepare_hdf5_inputs()
            self.sipros_search(config)
