import os
import sys
import configparser
import argparse
from argparse import Namespace


class SIPROSWorkflow:
    def __init__(self) -> None:
        # Default paths for third-party tools
        self.default_paths: dict[str, str] = {
            'raxport': 'tools/raxport',
            'sipros': 'tools/sipros',
            'feature_extractor': 'tools/aerith',
            'filter': 'tools/percolator',
            'deepfilter': 'tools/deepfilter',
            'assembly': 'tools/philosopher',
            'metaLP': 'tools/metaLP',
            'quantification': 'tools/ionquant'
        }

        # Initialize paths from config file if available
        self.cfg_file = 'workflow.cfg'
        self.paths: dict[str, str] = self.load_paths()

    def load_paths(self) -> dict[str, str]:
        paths: dict[str, str] = self.default_paths.copy()
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
                            help="Input raw/ft/mzml file path or directory, e.g., 'data/raw', 'A.raw,B.raw'")
        parser.add_argument('-e', '--element', required=True,
                            help="SIP label element, e.g., C, H, N, O, S")
        parser.add_argument('-f', '--fasta', required=True,
                            help="fasta file path")
        parser.add_argument('-r', '--range', required=True,
                            help="SIP label range, e.g., 0-100")
        parser.add_argument('-o', '--output', required=True,
                            help="Output directory path")

        args: Namespace = parser.parse_args()
        return args

    def run(self) -> None:
        args: Namespace = self.parse_arguments()

        # Call the feature extraction tool
        feature_cmd = f"{self.paths['feature_extractor']} {args.input} {args.label} {args.range} {args.output}"
        os.system(feature_cmd)

        # Call the filter tools
        filter_cmd = f"{self.paths['filter']} {args.input} {args.output}"
        os.system(filter_cmd)

        deepfilter_cmd = f"{self.paths['deepfilter']} {args.input} {args.output}"
        os.system(deepfilter_cmd)

        # Call the assembly tools
        assembly_cmd = f"{self.paths['assembly']} {args.output}"
        os.system(assembly_cmd)

        metaLP_cmd = f"{self.paths['metaLP']} {args.output}"
        os.system(metaLP_cmd)

        # Call the quantification tool
        quant_cmd = f"{self.paths['quantification']} {args.input} {args.output}"
        os.system(quant_cmd)

        # Output paths
        print(f"Results saved in {args.output}")


if __name__ == "__main__":
    workflow = SIPROSWorkflow()
    if len(sys.argv) == 1:
        print("SIPROS Workflow: A tool for integrating various parts of SIPROS into a complete workflow.")
        print("Use -h or --help to display help message.")
        sys.exit(0)
    workflow.run()
