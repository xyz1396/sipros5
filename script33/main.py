import os
import sys
import argparse
import csv
import logging
from logging import Logger
import time
from argparse import Namespace
from search import search
from filter import filter
from thread_allocation import (
    MIN_SIPROS_THREADS,
    available_cpu_count,
    effective_thread_count,
)
import warnings


class SIPROSWorkflow:
    def __init__(self) -> None:
        script_path = os.path.abspath(__file__)
        upper_path = os.path.dirname(os.path.dirname(script_path))
        self.upper_path = upper_path
        # Use the native tool set for the host platform. Keep all tools beside
        # their runtime DLLs and DIA-NN models in the repository/package tools
        # directory.
        if os.name == 'nt':
            raxport_name = 'Raxport-win-x64.exe'
            sipros_name = 'sipros.exe'
            aerith_name = 'aerith.exe'
        else:
            raxport_name = 'Raxport-linux-x64'
            sipros_name = 'sipros'
            aerith_name = 'aerith'
        self.defaultToolsPaths: dict[str, str] = {
            'raxport': os.path.join(upper_path, 'tools', raxport_name),
            'sipros': os.path.join(upper_path, 'tools', sipros_name),
            'aerith': os.path.join(upper_path, 'tools', aerith_name),
        }
        self.toolsPaths: dict[str, str] = self.defaultToolsPaths.copy()
        self.args: Namespace = self.parse_arguments()
        if not os.path.exists(self.args.output):
            os.makedirs(self.args.output)
        else:
            warnings.warn(message=f'{self.args.output} exists and will be overwritten')
        self.logger: Logger = self.initLogger(self.args.output)

    def parse_arguments(self) -> Namespace:
        epilog = """
Label search command demo : 
siproswf -i raw -f Ecoli.fasta -e C13 -t 40 -o wf_output
Regular search command demo : 
siproswf -i raw -f Ecoli.fasta -t 40 -o wf_output

citation:
1. Xiong, Y., Mueller, R.S., Feng, S., Guo, X. and Pan, C., 2024. Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins. Microbiome, 12.
2. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing. Bioinformatics Advances, vbae024
3. Guo, X., Li, Z., Yao, Q., Mueller, R.S., Eng, J.K., Tabb, D.L., Hervey IV, W.J. and Pan, C., 2018. Sipros ensemble improves database searching and filtering for complex metaproteomics. Bioinformatics, 34(5), pp.795-802
4. Wang, Y., Ahn, T.H., Li, Z. and Pan, C., 2013. Sipros/ProRata: a versatile informatics system for quantitative community proteomics. Bioinformatics, 29(16), pp.2064-2065
5. Pan, C., Kora, G., McDonald, W.H., Tabb, D.L., VerBerkmoes, N.C., Hurst, G.B., Pelletier, D.A., Samatova, N.F. and Hettich, R.L., 2006. ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation. Analytical chemistry, 78(20), pp.7121-7131
        """
        parser = argparse.ArgumentParser(
            description="sipros Workflow", prog="siproswf", epilog=epilog,
            formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i', '--input', required=True,
                            help="Input raw/.d/.d.zip/HDF5 file path or directory, e.g., 'data/raw', 'A.raw,B.h5,C.raw'")
        parser.add_argument(
            '--ptm', action='append', default=None,
            metavar='NAME_OR_TOKEN|default|none|all',
            help=("Replace the variable-PTM selection for FASTA search. Repeat "
                  "the option to select multiple PTMs; include 'default' to add "
                  "the profile defaults. 'none' disables all and 'all' selects "
                  "every compatible PTM. If omitted, profile defaults are used."),
        )
        parser.add_argument(
            '--fixed-ptm', action='append', default=None,
            metavar='NAME|default|none|all',
            help=("Replace the fixed-PTM selection for FASTA search or spectra-"
                  "library generation. Repeat to select multiple fixed PTMs; "
                  "'none' closes fixed CAM. If omitted, profile defaults are used."),
        )
        parser.add_argument(
            '--max-ptm-count', type=int, default=None, metavar='N',
            help=("Maximum variable PTMs allowed per peptide for FASTA search. "
                  "If omitted, the compiled search-profile default is used."),
        )
        parser.add_argument('--toleranceMS1', required=False, type=float, default=0.01,
            help="MS1 mass tolerance in Da (default: 0.01).")
        parser.add_argument('--toleranceMS2', required=False, type=float, default=0.01,
            help="MS2 mass tolerance in Da (default: 0.01). For ion trap MS2 data, 0.1 is recommended."
        )
        parser.add_argument('-e', '--element', required=False,
                            help="SIP isotope: C13, H2, N15, O18, or S34. Do not provide this flag for regular search")
        parser.add_argument('-r', '--range', required=False,
                            help=("SIP label percentage, range, or comma-separated "
                                  "list, e.g., 50, 49-51, or 1,2,3,49,50,51. "
                                  "Don't provide this flag for regular search"))
        parser.add_argument('-p', '--precision', required=False,
                            help="SIP label precision in percentage, e.g., 1. Don't provide this flag for regular search")
        parser.add_argument('--psm-tsv', required=False,
                            help="PSM TSV used to generate search-spectra SFI spectra libraries")
        parser.add_argument('--unlabeled-input', required=False,
                            help="Unlabeled raw/.d/.d.zip/HDF5 input used to generate search-spectra spectra libraries")
        parser.add_argument('--spectra-dir', required=False,
                            help="Reuse an existing SFI spectra-library directory instead of generating one")
        parser.add_argument(
            '--fast-sip-search', action='store_true',
            help=("Run regular target/decoy FASTA search and Aerith filtering, "
                  "build one target and one decoy SFI from the resulting "
                  "*_filtered_psms.tsv files, then run SIP spectra search"),
        )
        parser.add_argument('-f', '--fasta', required=True,
                            help="Fasta file path")
        parser.add_argument('-n', '--nPrecursor', required=False,
                            type=int, default=6,
                            help=("Raxport isotope-envelope apex m/z values selected per MSn scan "
                                  "before charge expansion (default: 6; consider 15 for DIA)"))
        parser.add_argument('--product-top-isotopes', required=False,
                            type=int, default=5,
                            help=("Theoretical isotope peaks retained for each "
                                  "predicted SIP product ion (default: 5)"))
        parser.add_argument('-t', '--thread', required=False, type=int, default=0,
                            help=("Total CPU-thread budget for the whole workflow "
                                  "(default: all CPUs available to this process). "
                                  "Regular FASTA target/decoy searches split the "
                                  "budget between two processes; other Sipros and "
                                  "Sipros search jobs receive at least 8 threads when "
                                  "the budget permits"))
        parser.add_argument(
            '--aerith-sample-parallelism', type=int, default=3,
            help=("Number of samples Aerith processes concurrently "
                  "(default: 3; higher values increase RAM use)"),
        )
        parser.add_argument('--topN', '--top-psms-per-scan', dest='topN', required=False, type=int, default=20,
                            help="Top PSM rows retained per scan for target and decoy searches before merge (default: 20)")
        parser.add_argument('--rt-tolerance', required=False, type=float, default=5.0,
                            help="Retention-time candidate window in minutes for spectra search (default: 5)")
        parser.add_argument('--sfi-envelope-top-n', required=False, type=int, default=3,
                            help=("Isotope peaks stored per precursor/product envelope "
                                  "and used by the compact SFI gate/MVH/XCorr stages "
                                  "(default: 3); final WDP regenerates full envelopes"))
        parser.add_argument('--mvh-cascade-top-n', required=False, type=int, default=150,
                            help=("Spectra-search MVH candidates retained per scan for "
                                  "XCorr/WDP scoring (default: 150)"))
        parser.add_argument('-o', '--output', required=True, help="Output directory path")
        parser.add_argument('--ignorePCT', action='store_true', 
                            help='Ignore SIP abundance features when filtering')
        parser.add_argument('--negative_control', required=False, type=str,
                            help="Negative control file name without extension, e.g., 'A' for 'A.raw', 'A,B' for 'A.raw,B.raw'.\n"
                                 "These files will be used to filter out false positive SIP-labeled PSMs")
        parser.add_argument('--label_threshold', required=False, type=float, default=2.0,
                            help="SIP label threshold in '%%' for filtering out false positive SIP-labeled PSMs")
        parser.add_argument('--dryrun', action='store_true', help='Run in dry run mode for test')

        args: Namespace = parser.parse_args()
        
        # Validate thread number
        if args.thread < 0:
            parser.error("Thread number must be non-negative (0 for all threads, or a positive integer)")
        if args.aerith_sample_parallelism <= 0:
            parser.error('--aerith-sample-parallelism must be a positive integer')
        available_threads = available_cpu_count()
        if args.thread > available_threads:
            warnings.warn(
                f"Requested {args.thread} threads, but only {available_threads} CPUs are available; "
                f"using {available_threads} threads"
            )
        args.thread = effective_thread_count(args.thread, available_threads)
        if args.thread < MIN_SIPROS_THREADS:
            warnings.warn(
                f"The {MIN_SIPROS_THREADS}-thread minimum for "
                f"Sipros cannot fit within a {args.thread}-thread "
                "workflow budget; those jobs will run serially using the full budget"
            )
        if args.topN <= 0:
            parser.error('--topN must be a positive integer')
        if args.nPrecursor <= 0:
            parser.error('--nPrecursor must be a positive integer')
        if args.product_top_isotopes <= 0:
            parser.error('--product-top-isotopes must be a positive integer')
        if args.rt_tolerance < 0:
            parser.error('--rt-tolerance must be non-negative')
        if args.sfi_envelope_top_n <= 0:
            parser.error('--sfi-envelope-top-n must be a positive integer')
        if args.mvh_cascade_top_n <= 0:
            parser.error('--mvh-cascade-top-n must be a positive integer')
        if args.max_ptm_count is not None and args.max_ptm_count < 0:
            parser.error('--max-ptm-count must be a non-negative integer')
        
        spectra_mode = bool(args.psm_tsv or args.unlabeled_input or args.spectra_dir)
        if args.fast_sip_search and spectra_mode:
            parser.error(
                "--fast-sip-search creates its own filtered PSMs and spectra "
                "library; do not combine it with --psm-tsv, "
                "--unlabeled-input, or --spectra-dir"
            )
        if spectra_mode and (
                args.ptm is not None or args.max_ptm_count is not None):
            parser.error(
                "--ptm and --max-ptm-count apply only to FASTA search; "
                "search-spectra uses the PTMs already encoded in its spectra library"
            )
        if args.spectra_dir and args.fixed_ptm is not None:
            parser.error(
                "--fixed-ptm cannot be used with --spectra-dir; the reused "
                "spectra library's chemistry metadata is authoritative"
            )
        if not args.element:
            args.element = "C13" if (spectra_mode or args.fast_sip_search) else "R"
        elif args.element != "R":
            normalized_element = args.element[0].upper() + args.element[1:]
            supported_isotopes = {"C13", "H2", "N15", "O18", "S34"}
            if normalized_element not in supported_isotopes:
                parser.error(
                    "--element must be one of C13, H2, N15, O18, or S34")
            args.element = normalized_element
        if spectra_mode and not args.spectra_dir and (not args.psm_tsv or not args.unlabeled_input):
            parser.error("search-spectra mode requires --psm-tsv and --unlabeled-input unless --spectra-dir is provided")
        if args.fast_sip_search and args.element == "R":
            parser.error("--fast-sip-search requires a SIP isotope such as C13")
        return args

    def initLogger(self, outputPath: str) -> Logger:
        logger: Logger = logging.getLogger('sipros_workflow')
        logger.setLevel(logging.INFO)
        file_handler = logging.FileHandler(
            f'{outputPath}/sipros_workflow.log', mode='w')
        stream_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)
        logger.info('sipros_workflow begin')
        return logger

    def make_search(self, *, element: str, output: str,
                    psm_tsv: str | None = None,
                    unlabeled_input: str | None = None,
                    spectra_dir: str | None = None,
                    stage_hdf5_copies: bool = True) -> search:
        return search(element=element,
                      toleranceMS1=self.args.toleranceMS1,
                      toleranceMS2=self.args.toleranceMS2,
                      sipRange=self.args.range,
                      step=self.args.precision,
                      raxportPath=self.toolsPaths['raxport'],
                      siprosPath=self.toolsPaths['sipros'],
                      fastaPath=self.args.fasta,
                      inputPath=self.args.input,
                      outputPath=output,
                      negative_control=self.args.negative_control,
                      threadNumber=int(self.args.thread),
                      logger=self.logger,
                      nPrecursor=self.args.nPrecursor,
                      dryrun=self.args.dryrun,
                      psmTsv=psm_tsv,
                      unlabeledInput=unlabeled_input,
                      spectraDir=spectra_dir,
                      topPsmsPerScan=self.args.topN,
                      ptms=self.args.ptm,
                      fixedPtms=self.args.fixed_ptm,
                      maxPtmCount=self.args.max_ptm_count,
                      rtTolerance=self.args.rt_tolerance,
                      sfiEnvelopeTopN=self.args.sfi_envelope_top_n,
                      mvhCascadeTopN=self.args.mvh_cascade_top_n,
                      stageHdf5Copies=stage_hdf5_copies)

    def make_filter(self, sipros_search: search, output: str, *,
                    assemble_proteins: bool, sip_isotope: str,
                    negative_control: str = "",
                    include_spectra: bool = False,
                    prediction_cache_path: str = "") -> filter:
        spectra_paths = None
        if include_spectra:
            spectra_paths = [
                sipros_search.hdf5_paths.get(
                    name, sipros_search.expected_hdf5_path(name)
                )
                for name in sipros_search.base_names
            ]
        return filter(baseNames=sipros_search.base_names,
                      outputPath=output,
                      aerithPath=self.toolsPaths['aerith'],
                      threadNumber=sipros_search.threadNumber,
                      logger=self.logger,
                      decoyPrefix=sipros_search.decoyPrefix,
                      ignorePCT=self.args.ignorePCT,
                      dryrun=self.args.dryrun,
                      fastaPath=self.args.fasta,
                      decoyPath=sipros_search.decoyPath,
                      assembleProteins=assemble_proteins,
                      element=sip_isotope,
                      negative_control=negative_control,
                      label_threshold=self.args.label_threshold,
                      fixedCam=not any(
                          value.strip().lower() == "none"
                          for value in (self.args.fixed_ptm or [])
                      ),
                      sipIsotope=sip_isotope,
                      ptms=self.args.ptm,
                      fixedPtms=self.args.fixed_ptm,
                      maxPtmCount=self.args.max_ptm_count,
                      productTopIsotopes=self.args.product_top_isotopes,
                      quantTopIsotopes=self.args.nPrecursor,
                      sampleParallelism=self.args.aerith_sample_parallelism,
                      spectraPaths=spectra_paths,
                      predictionCachePath=prediction_cache_path)

    def report_filtered_psms(self, regular_search: search,
                             regular_output: str) -> None:
        for base_name in regular_search.base_names:
            path = os.path.join(
                regular_output, base_name, f'{base_name}_filtered_psms.tsv'
            )
            accepted = 0
            if os.path.exists(path):
                with open(path, newline='') as stream:
                    reader = csv.DictReader(stream, delimiter='\t')
                    for row in reader:
                        label = (row.get('Label') or '1').strip()
                        if label == '1':
                            accepted += 1
            if accepted == 0:
                self.logger.warning(
                    f'Aerith regular-search filtering retained 0 target PSMs '
                    f'for {base_name}; this sample contributes no target '
                    'spectra to the SFI library'
                )
            else:
                self.logger.info(
                    f'Aerith regular-search filtering retained {accepted:,} '
                    f'target PSMs for {base_name}'
                )

    def run_fast_sip_search(self) -> search:
        regular_output = os.path.join(self.args.output, 'regular')
        spectra_output = os.path.join(self.args.output, 'spectra_search')
        prediction_cache = os.path.join(
            regular_output, 'regular_search_predictions'
        )
        for path in (regular_output, spectra_output):
            os.makedirs(path, exist_ok=True)

        phase_timings: list[tuple[str, float]] = []
        phase_started = time.perf_counter()
        self.logger.info('Fast SIP phase 1/4: regular target/decoy FASTA search')
        regular_search = self.make_search(
            element='R', output=regular_output, stage_hdf5_copies=False
        )
        regular_search.run()
        phase_timings.append(
            ('1/4 Regular FASTA search', time.perf_counter() - phase_started)
        )

        phase_started = time.perf_counter()
        self.logger.info('Fast SIP phase 2/4: Aerith filtering of regular PSMs')
        if not self.args.dryrun:
            for suffix in ('.bin', '.spectrum', '.rt'):
                cache_file = prediction_cache + suffix
                if os.path.exists(cache_file):
                    os.remove(cache_file)
        self.logger.info(
            'Aerith regular filtering predicts DIA-NN spectra and RT once and '
            'caches every unique target PIN peptide-charge form in the same run'
        )
        regular_filter = self.make_filter(
            regular_search, regular_output,
            assemble_proteins=True, sip_isotope='',
            include_spectra=True,
            prediction_cache_path=prediction_cache,
        )
        regular_filter.run()
        if not self.args.dryrun:
            self.report_filtered_psms(regular_search, regular_output)
        phase_timings.append(
            ('2/4 Regular Aerith filter', time.perf_counter() - phase_started)
        )

        phase_started = time.perf_counter()
        self.logger.info(
            'Fast SIP phase 3/4: filtered-PSM SFI generation and spectra search'
        )
        spectra_search = self.make_search(
            element=self.args.element,
            output=spectra_output,
            psm_tsv=regular_output,
            unlabeled_input=regular_output,
            stage_hdf5_copies=False,
        )
        spectra_search.base_names = list(regular_search.base_names)
        spectra_search.hdf5_paths = dict(regular_search.hdf5_paths)
        spectra_search.generatedSpectraDir = spectra_output
        spectra_search.reverse_fasta_sequences()
        generated = spectra_output
        if not self.args.dryrun:
            generated = spectra_search.generate_or_reuse_spectra_library()
            spectra_search.search_spectra_samples(generated)
        phase_timings.append(
            ('3/4 SFI + spectra search', time.perf_counter() - phase_started)
        )

        phase_started = time.perf_counter()
        self.logger.info('Fast SIP phase 4/4: Aerith filtering and reporting')
        spectra_filter = self.make_filter(
            spectra_search, spectra_output,
            assemble_proteins=True,
            sip_isotope=self.args.element,
            negative_control=self.args.negative_control or '',
            include_spectra=True,
            prediction_cache_path=prediction_cache,
        )
        self.logger.info(
            'Aerith spectra-search filtering reuses cached target predictions; '
            'generated-decoy spectra and RT are predicted in memory and are '
            'not written to the target-only cache'
        )
        spectra_filter.run()
        phase_timings.append(
            ('4/4 SIP Aerith + reports', time.perf_counter() - phase_started)
        )
        width = max(len(label) for label, _ in phase_timings)
        timing_lines = ['Fast SIP phase timing (wall clock)']
        for label, seconds in phase_timings:
            timing_lines.append(f'  {label:<{width}} : {seconds:9.3f} s')
        timing_lines.append(
            f'  {"Fast SIP total":<{width}} : '
            f'{sum(seconds for _, seconds in phase_timings):9.3f} s'
        )
        self.logger.info('\n'.join(timing_lines))
        return spectra_search

    def run(self) -> None:
        start_time = time.perf_counter()
        spectra_mode = bool(self.args.psm_tsv or self.args.unlabeled_input or self.args.spectra_dir)
        if self.args.fast_sip_search:
            self.run_fast_sip_search()
            running_time = time.perf_counter() - start_time
            self.logger.info(f'All job done. Results are in {self.args.output}.')
            self.logger.info(f'Total wall time: {running_time:.3f} s')
            return

        sipros_search = self.make_search(
            element=self.args.element,
            output=self.args.output,
            psm_tsv=self.args.psm_tsv,
            unlabeled_input=self.args.unlabeled_input,
            spectra_dir=self.args.spectra_dir,
        )
        if spectra_mode:
            sipros_search.run_search_spectra()
        else:
            sipros_search.run()

        sipros_filter = self.make_filter(
            sipros_search, self.args.output,
            assemble_proteins=not spectra_mode,
            sip_isotope=self.args.element or '',
            negative_control=self.args.negative_control or '',
            include_spectra=not spectra_mode,
        )
        sipros_filter.run()

        running_time = time.perf_counter() - start_time
        self.logger.info(f'All job done. Results are in {self.args.output}.')
        self.logger.info(f'Total wall time: {running_time:.3f} s')


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("SIPROS Workflow: A tool for integrating various parts of SIPROS into a complete workflow.")
        print("Use -h or --help to display help message.")
        sys.exit(0)
    workflow = SIPROSWorkflow()
    workflow.run()
