import os
import sys
import configparser
import argparse
import logging
from logging import Logger
import time
from argparse import Namespace
from search import search
import warnings

class SIPROSWorkflow:
    def __init__(self) -> None:
        script_path = os.path.abspath(__file__)
        upper_path = os.path.dirname(os.path.dirname(script_path))
        # Default paths for tools
        self.defaultToolsPaths: dict[str, str] = {
            'configTemplates': f'{upper_path}/configTemplates',
            'configGenerator': f'{upper_path}/tools/configGenerator',
            'raxport': f'{upper_path}/tools/raxport',
            'sipros': f'{upper_path}/tools/sipros',
            'feature_extractor': f'{upper_path}/tools/aerith',
            'filter': f'{upper_path}/tools/percolator',
            'deepfilter': f'{upper_path}/tools/deepfilter',
            'assembly': f'{upper_path}/tools/philosopher',
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
            warnings.warn(f'{self.args.output} exists and will be overwritten')
        self.logger: Logger = self.initLogger(self.args.output)

    def load_paths(self) -> dict[str, str]:
        paths: dict[str, str] = self.defaultToolsPaths.copy()
        if os.path.exists(path=self.cfg_file):
            config = configparser.ConfigParser()
            config.read(filenames=self.cfg_file)
            for key in paths:
                if config.has_option(section='Paths', option=key):
                    paths[key] = config.get(section='Paths', option=key)
        return paths

    def parse_arguments(self) -> Namespace:
        citation = """
citation:
1. Xiong, Y., Mueller, R.S., Feng, S., Guo, X. and Pan, C., 2024. Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins. Microbiome, 12.
2. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing. Bioinformatics Advances, vbae024
3. Guo, X., Li, Z., Yao, Q., Mueller, R.S., Eng, J.K., Tabb, D.L., Hervey IV, W.J. and Pan, C., 2018. Sipros ensemble improves database searching and filtering for complex metaproteomics. Bioinformatics, 34(5), pp.795-802
4. Wang, Y., Ahn, T.H., Li, Z. and Pan, C., 2013. Sipros/ProRata: a versatile informatics system for quantitative community proteomics. Bioinformatics, 29(16), pp.2064-2065
5. Pan, C., Kora, G., McDonald, W.H., Tabb, D.L., VerBerkmoes, N.C., Hurst, G.B., Pelletier, D.A., Samatova, N.F. and Hettich, R.L., 2006. ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation. Analytical chemistry, 78(20), pp.7121-7131
        """
        parser = argparse.ArgumentParser(
            description="sipros Workflow", prog="siproswf", epilog=citation,
            formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i', '--input', required=True,
                            help="Input raw/mzml file path or directory, e.g., 'data/raw', 'A.raw,B.raw'")
        parser.add_argument('-e', '--element', required=False,
                            help="SIP label element, e.g., C13, H2, N15, O18, S33, S34")
        parser.add_argument('-r', '--range', required=False,
                            help="SIP label range, e.g., 0-100")
        parser.add_argument('-p', '--precision', required=False,
                            help="SIP label precision in percentage, e.g., 1")
        parser.add_argument('-f', '--fasta', required=True,
                            help="fasta file path")
        parser.add_argument('-s', '--split_FT2_file', required=False,
                            type=int, nargs='?', const=20000,
                            help="scans number in splitted FT2 file, no split in default")
        parser.add_argument('-t', '--thread', required=False,
                            help="thread number to be limited, all threads in default")
        parser.add_argument('-o', '--output', required=True,
                            help="Output directory path")

        args: Namespace = parser.parse_args()
        return args
    
    def initLogger(self, outputPath: str) -> Logger:
        logger: Logger = logging.getLogger('sipros_workflow')
        logger.setLevel(logging.INFO)
        file_handler = logging.FileHandler(f'{outputPath}/sipros_workflow.log', mode='w')
        stream_handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)
        logger.info('sipros_workflow begin')
        return logger

    def run(self) -> None:
        start_time: float = time.time()
        
        # run SIPROS search
        sipros_search = search(element=self.args.element, 
                               sipRange=self.args.range,
                               step=self.args.precision,
                               configTemplatePath=self.toolsPaths['configTemplates'],
                               configGeneratorPath=self.toolsPaths['configGenerator'],
                               raxportPath=self.toolsPaths['raxport'],
                               scansPerFT2=self.args.split_FT2_file, 
                               siprosPath=self.toolsPaths['sipros'],
                               fastaPath=self.args.fasta, 
                               inputPath=self.args.input, 
                               outputPath=self.args.output,
                               threadNumber=int(self.args.thread), 
                               logger=self.logger)
        sipros_search.run()
        
        # Call the feature extraction tool


        # Call the filter tools

        # Call the assembly tools

        # Call the quantification tool

        # Output paths
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
