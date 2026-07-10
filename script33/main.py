import os
import sys
import configparser
import argparse
import logging
from logging import Logger
import time
from argparse import Namespace
from search import search
from filter import filter
from assembly import assembly
import warnings


class SIPROSWorkflow:
    def __init__(self) -> None:
        script_path = os.path.abspath(__file__)
        upper_path = os.path.dirname(os.path.dirname(script_path))
        self.upper_path = upper_path
        # Default paths for tools
        self.defaultToolsPaths: dict[str, str] = {
            'configTemplates': f'{upper_path}/configTemplates',
            'raxport': f'{upper_path}/tools/Raxport-linux-x64',
            'sipros': f'{upper_path}/tools/sipros',
            'filter': f'{upper_path}/tools/percolator',
            'deepfilter': f'{upper_path}/tools/deepfilter',
            'assembly': f'{upper_path}/tools/philosopher-v5.1.2',
            'metaLP': f'{upper_path}/tools/metaLP',
            'quantification': f'{upper_path}/tools/ionquant'
        }

        # Initialize paths from config file if available
        self.cfg_file = 'workflow.cfg'
        self.toolsPaths: dict[str, str] = self.load_paths()
        self.args: Namespace = self.parse_arguments()
        if not os.path.exists(self.args.output):
            os.makedirs(self.args.output)
        else:
            warnings.warn(message=f'{self.args.output} exists and will be overwritten')
        self.logger: Logger = self.initLogger(self.args.output)

    def load_paths(self) -> dict[str, str]:
        paths: dict[str, str] = self.defaultToolsPaths.copy()
        cfg_path = os.path.join(self.upper_path, self.cfg_file)
        if os.path.exists(path=cfg_path):
            config = configparser.ConfigParser()
            config.read(filenames=cfg_path)
            for key in paths:
                if config.has_option(section='Paths', option=key):
                    paths[key] = config.get(section='Paths', option=key)
        return paths

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
        parser.add_argument('--sPTM', required=False, type=str, help=(
            "Static PTM for regular search. Example: 'C|5,8,2,2,0,1'\n"
            "This means a static modification (e.g., alkylation) on residue C, "
            "with new element composition: C,H,O,N,P,S = 5,8,2,2,0,1."
        ))
        parser.add_argument('--vPTM', required=False, type=str, help=(
            "Variable PTM for regular search, e.g., '~|M|0,0,1,0,0,0;!|NQ|0,-1,1,-1,0,0'.\n"
            "This means a variable modification on residue M (Oxidation or Hydroxylation, element composition change: 0,0,1,0,0,0 for C,H,O,N,P,S),\n"
            "and on N/Q (Deamidation or Citrullination, element composition change: 0,-1,1,-1,0,0 for C,H,O,N,P,S)."
        ))
        parser.add_argument('--toleranceMS1', required=False, type=float, default=0.01,
            help="MS1 mass tolerance in Da (default: 0.01).")
        parser.add_argument('--toleranceMS2', required=False, type=float, default=0.01,
            help="MS2 mass tolerance in Da (default: 0.01). For ion trap MS2 data, 0.1 is recommended."
        )
        parser.add_argument('-e', '--element', required=False,
                            help="SIP isotope: C13, H2, N15, O18, or S34. Do not provide this flag for regular search")
        parser.add_argument('-r', '--range', required=False,
                            help="SIP label range, e.g., 0-100. Don't provide this flag for regular search")
        parser.add_argument('-p', '--precision', required=False,
                            help="SIP label precision in percentage, e.g., 1. Don't provide this flag for regular search")
        parser.add_argument('--psm-tsv', required=False,
                            help="PSM TSV used to generate search-spectra HDF5 spectra libraries")
        parser.add_argument('--unlabeled-input', required=False,
                            help="Unlabeled raw/.d/.d.zip/HDF5 input used to generate search-spectra spectra libraries")
        parser.add_argument('--spectra-dir', required=False,
                            help="Reuse an existing HDF5 spectra-library directory instead of generating one")
        parser.add_argument('-f', '--fasta', required=True,
                            help="Fasta file path")
        parser.add_argument('-n', '--nPrecursor', required=False,
                            type=int, default=6,
                            help=("Raxport isotope-envelope apex m/z values selected per MSn scan "
                                  "before charge expansion (default: 6; consider 15 for DIA)"))
        parser.add_argument('-t', '--thread', required=False, type=int, default=0,
                            help="Thread number to be limited, all threads in default")
        parser.add_argument('--topN', '--top-psms-per-scan', dest='topN', required=False, type=int, default=8,
                            help="Top PSM rows retained per scan for target and decoy searches before merge (default: 8)")
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
        if args.thread == 0:
            args.thread = os.cpu_count() or 1  # Use all available CPUs
        if args.topN <= 0:
            parser.error('--topN must be a positive integer')
        if args.nPrecursor <= 0:
            parser.error('--nPrecursor must be a positive integer')
        
        spectra_mode = bool(args.psm_tsv or args.unlabeled_input or args.spectra_dir)
        if not args.element:
            args.element = "C13" if spectra_mode else "R"
        elif args.element != "R":
            normalized_element = args.element[0].upper() + args.element[1:]
            supported_isotopes = {"C13", "H2", "N15", "O18", "S34"}
            if normalized_element not in supported_isotopes:
                parser.error(
                    "--element must be one of C13, H2, N15, O18, or S34")
            args.element = normalized_element
        if spectra_mode and not args.spectra_dir and (not args.psm_tsv or not args.unlabeled_input):
            parser.error("search-spectra mode requires --psm-tsv and --unlabeled-input unless --spectra-dir is provided")
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

    def run(self) -> None:
        start_time: float = time.time()
        spectra_mode = bool(self.args.psm_tsv or self.args.unlabeled_input or self.args.spectra_dir)

        sipros_search = search(element=self.args.element,
                               toleranceMS1=self.args.toleranceMS1,
                               toleranceMS2=self.args.toleranceMS2,
                               sipRange=self.args.range,
                               step=self.args.precision,
                               configTemplatePath=self.toolsPaths['configTemplates'],
                               raxportPath=self.toolsPaths['raxport'],
                               siprosPath=self.toolsPaths['sipros'],
                               fastaPath=self.args.fasta,
                               inputPath=self.args.input,
                               outputPath=self.args.output,
                               negative_control=self.args.negative_control,
                               threadNumber=int(self.args.thread),
                               logger=self.logger,
                               nPrecursor=self.args.nPrecursor,
                               dryrun=self.args.dryrun,
                               psmTsv=self.args.psm_tsv,
                               unlabeledInput=self.args.unlabeled_input,
                               spectraDir=self.args.spectra_dir,
                               topPsmsPerScan=self.args.topN)
        if spectra_mode:
            sipros_search.run_search_spectra()
        else:
            sipros_search.run()

        sipros_filter = filter(baseNames=sipros_search.base_names,
                               outputPath=self.args.output,
                               percolatorPath=self.toolsPaths['filter'],
                               threadNumber=sipros_search.threadNumber,
                               logger=self.logger,
                               ignorePCT=self.args.ignorePCT,
                               dryrun=self.args.dryrun)
        sipros_filter.run()

        sipros_assembly = assembly(baseNames=sipros_search.base_names,
                                   philosopherPath=self.toolsPaths['assembly'],
                                   percolatorPath=self.toolsPaths['filter'],
                                   fastaPath=self.args.fasta,
                                   decoyPath=sipros_search.decoyPath,
                                   outputPath=self.args.output,
                                   threadNumber=sipros_search.threadNumber,
                                   negative_control=self.args.negative_control,
                                   element=self.args.element,
                                   label_threshold=self.args.label_threshold,
                                   logger=self.logger,
                                   decoyPrefix=sipros_search.decoyPrefix)
        sipros_assembly.run()

        end_time = time.time()
        running_time = end_time - start_time
        self.logger.info(f'All job done. Results are in {self.args.output}.')
        self.logger.info(f'Total running time: {running_time} seconds')


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("SIPROS Workflow: A tool for integrating various parts of SIPROS into a complete workflow.")
        print("Use -h or --help to display help message.")
        sys.exit(0)
    workflow = SIPROSWorkflow()
    workflow.run()
