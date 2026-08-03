from logging import Logger
import os
import concurrent.futures
import shutil
import shlex
from pathlib import Path
from command_runner import run_logged_command
from thread_allocation import (
    MIN_SIPROS_THREADS,
    ThreadAllocation,
    allocate_threads,
    available_cpu_count,
    effective_thread_count,
    thread_env_updates,
)


class search:
    def __init__(self, toleranceMS1: float, toleranceMS2: float,
                 sipRange: str, step: str,
                 raxportPath: str, siprosPath: str, fastaPath: str,
                 inputPath: str, outputPath: str, negative_control: str,
                 threadNumber: int, logger: Logger, element="R",
                 nPrecursor=6, dryrun=False, psmTsv: str | None = None,
                 unlabeledInput: str | None = None, spectraDir: str | None = None,
                 topPsmsPerScan: int = 20, ptms: list[str] | None = None,
                 fixedPtms: list[str] | None = None,
                 maxPtmCount: int | None = None,
                 rtTolerance: float = 5.0,
                 sfiEnvelopeTopN: int = 3,
                 mvhCascadeTopN: int = 150,
                 stageHdf5Copies: bool = True) -> None:
        self.core_count = available_cpu_count()
        self.element = element
        self.toleranceMS1 = toleranceMS1
        self.toleranceMS2 = toleranceMS2
        self.sipRange = sipRange
        self.step = step
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
        self.ptms = None if ptms is None else list(ptms)
        self.fixedPtms = None if fixedPtms is None else list(fixedPtms)
        self.maxPtmCount = maxPtmCount
        self.rtTolerance = rtTolerance
        self.sfiEnvelopeTopN = sfiEnvelopeTopN
        self.mvhCascadeTopN = mvhCascadeTopN
        self.stageHdf5Copies = stageHdf5Copies
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
            def decoy_sequence(value: str) -> str:
                if self.element == 'R' and value:
                    # Regular search explicitly models intact and excised
                    # protein N termini.  Preserve the first target residue in
                    # the decoy so those terminal proteoforms have a symmetric
                    # target/decoy search space; SIP keeps its legacy reversal.
                    return value[0] + value[:0:-1]
                return value[::-1]
            for line in fasta:
                if line.startswith('>'):
                    if header:
                        output.write(header + '\n' + decoy_sequence(sequence) + '\n')
                    header = '>' + self.decoyPrefix + line[1:].strip()
                    sequence = ''
                else:
                    sequence += line.strip()
            if header:
                output.write(header + '\n' + decoy_sequence(sequence) + '\n')

    def input_entries(self, input_path: str) -> list[str]:
        if ',' in input_path:
            self.logger.info(f'{input_path} is a file list')
            return [
                item.strip()
                for item in input_path.split(',')
                if item.strip()
            ]
        path = Path(input_path)
        if self.is_raw_input(input_path) or self.is_hdf5_input(input_path):
            return [input_path]
        if path.is_dir():
            self.logger.info(f'{input_path} is a directory')
            return [str(path / name) for name in sorted(os.listdir(path))]
        return [input_path]

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
            if not self.stageHdf5Copies:
                if not os.path.lexists(expected) and not self.dryrun:
                    os.symlink(os.path.realpath(hdf5_file), expected)
                self.hdf5_paths[base] = expected
                self.logger.info(
                    f'Reusing HDF5 input through lightweight link: {expected}'
                )
                continue
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

    def direct_sip_args(self) -> str:
        if self.element == 'R':
            return ''
        sip_range = self.sipRange if self.sipRange is not None else '0-100'
        sip_step = self.step if self.step is not None else '1'
        return (f' -a {self.q(self.element)} -b {self.q(sip_range)}'
                f' -s {self.q(sip_step)}')

    def direct_ptm_args(self) -> str:
        args = ''.join(f' --ptm {self.q(ptm)}' for ptm in (self.ptms or []))
        args += self.fixed_ptm_args()
        if self.maxPtmCount is not None:
            args += f' --max-ptm-count {self.maxPtmCount}'
        return args

    def fixed_ptm_args(self) -> str:
        return ''.join(
            f' --fixed-ptm {self.q(ptm)}' for ptm in (self.fixedPtms or [])
        )

    def regular_fasta_search(self):
        """Prepare caches, then search each H5 with one target/decoy pair."""
        ptm_args = self.direct_ptm_args()
        tolerance_args = (
            f' --tolerance-ms1 {self.q(self.toleranceMS1)}'
            f' --tolerance-ms2 {self.q(self.toleranceMS2)}'
        )
        output_root = Path(self.outPutPath)
        output_root.mkdir(parents=True, exist_ok=True)

        chemistry_args = (
            f'{tolerance_args}{ptm_args}'
            ' --precursor-source ms1-neighborhood'
        )
        target_cache = output_root / 'target.sfi'
        decoy_cache = output_root / 'decoy.sfi'

        # Build/load both persistent caches first.  Unlike search jobs, these
        # two preparation processes may receive fewer than eight cores when
        # -t is small so they can remain concurrent while using the full
        # requested CPU budget.
        prepare_commands = [
            (
                f'{self.q(self.siprosPath)} search-fasta '
                f'-fasta {self.q(self.fastaPath)} --prepare-only'
                f'{chemistry_args} --fragment-index-cache '
                f'{self.q(target_cache)}'
            ),
            (
                f'{self.q(self.siprosPath)} search-fasta '
                f'-fasta {self.q(self.decoyPath)} --prepare-only'
                f'{chemistry_args} --fragment-index-cache '
                f'{self.q(decoy_cache)}'
            ),
        ]
        prepare_allocation = allocate_threads(
            self.threadNumber,
            len(prepare_commands),
        )
        self.log_thread_allocation(
            'Sipros Regular FASTA cache preparation', prepare_allocation
        )
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=prepare_allocation.worker_count) as executor:
            futures = [
                executor.submit(self.run_command_sipros, cmd, threads)
                for cmd, threads in zip(
                    prepare_commands, prepare_allocation.task_threads
                )
            ]
            for future in concurrent.futures.as_completed(futures):
                future.result()

        common_search_args = (
            f' --top-psms-per-scan {self.topPsmsPerScan}'
            f'{chemistry_args}'
        )
        search_pairs: list[tuple[str, str]] = []
        for base_name in self.base_names:
            scan_arg = f' -f {self.q(self.hdf5_paths[base_name])}'
            sample_dir = output_root / base_name
            sample_dir.mkdir(parents=True, exist_ok=True)
            target_pin_name = f'{base_name}_target.pin'
            decoy_pin_name = f'{base_name}_decoy.pin'
            search_pairs.append((
                (
                    f'{self.q(self.siprosPath)} search-fasta '
                    f'-fasta {self.q(self.fastaPath)}{scan_arg} '
                    f'-o {self.q(sample_dir)} '
                    f'--pin-output {self.q(target_pin_name)} --pin-label 1'
                    f'{common_search_args} --fragment-index-cache '
                    f'{self.q(target_cache)}'
                ),
                (
                    f'{self.q(self.siprosPath)} search-fasta '
                    f'-fasta {self.q(self.decoyPath)}{scan_arg} '
                    f'-o {self.q(sample_dir)} '
                    f'--pin-output {self.q(decoy_pin_name)} --pin-label -1'
                    f'{common_search_args} --fragment-index-cache '
                    f'{self.q(decoy_cache)}'
                ),
            ))

        # Regular search is most efficient with two cache-query processes.
        # Split the complete -t budget between the target and decoy for one
        # sample, wait for both, then advance to the next sample.  With an odd
        # budget the target receives the extra thread.  A one-thread budget
        # necessarily runs the pair sequentially rather than oversubscribing.
        search_allocation = allocate_threads(self.threadNumber, 2)
        self.log_thread_allocation(
            'Sipros Regular FASTA paired target/decoy cache-H5 search',
            search_allocation,
        )
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=search_allocation.worker_count) as executor:
            for commands in search_pairs:
                futures = [
                    executor.submit(
                        self.run_command_sipros,
                        command,
                        threads,
                    )
                    for command, threads in zip(
                        commands, search_allocation.task_threads
                    )
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()

    def sip_fasta_search(self):
        """Keep SIP FASTA assignment on Raxport's precursor candidates."""
        sip_args = self.direct_sip_args()
        ptm_args = self.direct_ptm_args()
        tolerance_args = (
            f' --tolerance-ms1 {self.q(self.toleranceMS1)}'
            f' --tolerance-ms2 {self.q(self.toleranceMS2)}'
        )
        commands: list[str] = []
        for base_name in self.base_names:
            hdf5_path = self.hdf5_paths[base_name]
            sample_dir = f'{self.outPutPath}/{base_name}'
            target_pin_name = f'{base_name}_target.pin'
            decoy_pin_name = f'{base_name}_decoy.pin'
            commands.append(
                f'{self.q(self.siprosPath)} search-fasta -fasta {self.q(self.fastaPath)} '
                f'-f {self.q(hdf5_path)} -o {self.q(sample_dir)} --pin-output {self.q(target_pin_name)}{sip_args} '
                f'--pin-label 1 --precursor-source raxport-candidates '
                f'--top-psms-per-scan {self.topPsmsPerScan}{tolerance_args}{ptm_args}'
            )
            commands.append(
                f'{self.q(self.siprosPath)} search-fasta -fasta {self.q(self.decoyPath)} '
                f'-f {self.q(hdf5_path)} -o {self.q(sample_dir)} --pin-output {self.q(decoy_pin_name)}{sip_args} '
                f'--pin-label -1 --precursor-source raxport-candidates '
                f'--top-psms-per-scan {self.topPsmsPerScan}{tolerance_args}{ptm_args}'
            )
        allocation = allocate_threads(
            self.threadNumber,
            len(commands),
            minimum_threads_per_task=MIN_SIPROS_THREADS,
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

    def sipros_search(self):
        if self.element == 'R':
            self.regular_fasta_search()
        else:
            self.sip_fasta_search()

    def resolve_or_convert_unlabeled_hdf5(self) -> str:
        if self.unlabeledInput is None or self.unlabeledInput == '':
            self.logger.error('--unlabeled-input is required when --spectra-dir is not provided')
            raise SystemExit(1)
        unlabeled_path = Path(self.unlabeledInput)
        if unlabeled_path.is_dir():
            real_directory = os.path.realpath(unlabeled_path)
            hdf5_files = [
                path for path in Path(real_directory).rglob('*')
                if path.is_file() and self.is_hdf5_input(str(path))
            ]
            if not hdf5_files:
                self.logger.error(
                    f'No HDF5 files found under unlabeled input directory: '
                    f'{real_directory}'
                )
                raise SystemExit(1)
            self.logger.info(
                f'Using {len(hdf5_files)} HDF5 files under {real_directory} '
                'to match filtered regular-search PSMs'
            )
            return real_directory
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

    def generate_or_reuse_spectra_library(self) -> str:
        def require_target_decoy_pair(directory: str) -> list[Path]:
            files = list(Path(directory).glob('*.sfi'))
            decoys = [path for path in files if 'decoy' in path.name.lower()]
            targets = [path for path in files if 'decoy' not in path.name.lower()]
            if len(targets) != 1 or len(decoys) != 1:
                self.logger.error(
                    f'SFI spectra library must contain exactly one target and '
                    f'one generated-decoy index; found {len(targets)} target '
                    f'and {len(decoys)} decoy files in {directory}'
                )
                raise SystemExit(1)
            return files

        if self.spectraDir:
            if self.fixedPtms:
                self.logger.error(
                    '--fixed-ptm cannot be used with --spectra-dir; the reused '
                    "spectra library's chemistry metadata is authoritative"
                )
                raise SystemExit(1)
            spectra_dir = os.path.realpath(self.spectraDir)
            if not os.path.isdir(spectra_dir):
                self.logger.error(f'--spectra-dir does not exist or is not a directory: {spectra_dir}')
                raise SystemExit(1)
            require_target_decoy_pair(spectra_dir)
            self.logger.info(f'Reusing SFI spectra library from {spectra_dir}')
            return spectra_dir
        if not self.psmTsv:
            self.logger.error('--psm-tsv is required to generate a spectra library')
            raise SystemExit(1)
        unlabeled_hdf5 = self.resolve_or_convert_unlabeled_hdf5()
        os.makedirs(self.generatedSpectraDir, exist_ok=True)
        sip_range = self.sipRange if self.sipRange is not None else '0-100'
        sip_step = self.step if self.step is not None else '1'
        cmd = (f'{self.q(self.siprosPath)} experimental-spectra '
               f'-i {self.q(self.psmTsv)} -f {self.q(unlabeled_hdf5)} '
               f'-o {self.q(self.generatedSpectraDir)} -a {self.element} '
               f'-b {self.q(sip_range)} -s {self.q(sip_step)} '
               f'--decoy -t {self.threadNumber}'
               f' --envelope-top-n {getattr(self, "sfiEnvelopeTopN", 3)}'
               f'{self.fixed_ptm_args()}')
        if not self.dryrun:
            self.run_command(cmd, threads=self.threadNumber)
            require_target_decoy_pair(self.generatedSpectraDir)
            non_sfi = [path for path in Path(self.generatedSpectraDir).iterdir()
                       if path.is_file() and path.suffix.lower() != '.sfi']
            if non_sfi:
                self.logger.error(f'Unexpected non-SFI spectra intermediates generated: {non_sfi}')
                raise SystemExit(1)
        return self.generatedSpectraDir

    def search_spectra_samples(self, spectra_dir: str):
        search_pairs: list[tuple[str, str]] = []
        for base_name in self.base_names:
            hdf5_path = self.hdf5_paths[base_name]
            sample_dir = f'{self.outPutPath}/{base_name}'
            os.makedirs(sample_dir, exist_ok=True)
            common = (
                f'{self.q(self.siprosPath)} search-spectra '
                f'-f {self.q(hdf5_path)} --sfi {self.q(spectra_dir)} '
                f'-o {self.q(sample_dir)} '
                f'--tolerance-ms1 {self.toleranceMS1} --tolerance-ms1-unit da '
                f'--tolerance-ms2 {self.toleranceMS2} --tolerance-ms2-unit da '
                f'--rt-tolerance {getattr(self, "rtTolerance", 5.0)} '
                f'--score-envelope-top-n {getattr(self, "sfiEnvelopeTopN", 3)} '
                f'--mvh-cascade-top-n {getattr(self, "mvhCascadeTopN", 150)} '
                f'--top-psms-per-scan {self.topPsmsPerScan}'
            )
            search_pairs.append((
                f'{common} --sfi-label target',
                f'{common} --sfi-label decoy',
            ))
        allocation = allocate_threads(self.threadNumber, 2)
        self.log_thread_allocation(
            'Sipros spectra paired target/decoy SFI-H5 search', allocation
        )
        if allocation.worker_count == 0:
            return
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=allocation.worker_count) as executor:
            for commands in search_pairs:
                futures = [
                    executor.submit(
                        self.run_command_sipros,
                        f'{command} -t {threads}',
                        threads,
                    )
                    for command, threads in zip(
                        commands, allocation.task_threads
                    )
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()

    def run_search_spectra(self) -> None:
        if self.element == 'R':
            self.logger.error('search-spectra mode requires a SIP element such as C13')
            raise SystemExit(1)
        self.logger.info(
            f'Workflow CPU allocation: {self.threadNumber} cores '
            f'({self.core_count} available)'
        )
        self.getInputFiles()
        self.create_sample_directories()
        if not self.dryrun:
            self.prepare_hdf5_inputs()
            spectra_dir = self.generate_or_reuse_spectra_library()
            self.search_spectra_samples(spectra_dir)

    def run(self) -> None:
        self.reverse_fasta_sequences()
        self.logger.info(
            f'Workflow CPU allocation: {self.threadNumber} cores '
            f'({self.core_count} available)'
        )
        self.getInputFiles()
        self.create_sample_directories()
        if not self.dryrun:
            self.prepare_hdf5_inputs()
            self.sipros_search()
