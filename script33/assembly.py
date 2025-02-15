from logging import Logger
import multiprocessing
import os
import subprocess
import concurrent.futures
import pandas as pd
from lxml import etree
import re
import shutil

class assembly:
    def __init__(self, baseNames: list[str], philosopherPath:str, fastaPath:str, decoyPath:str, outputPath: str,
                threadNumber: int, logger: Logger):
        self.philosopherPath = philosopherPath
        self.fastaPath = fastaPath
        self.decoyPath = decoyPath
        self.baseNames = baseNames
        self.outputPath = outputPath
        # merged fasta file for philosopher
        self.targetDecoyPath = f'{outputPath}/targetDecoy.faa'
        self.logger = logger
        self.threadNumber = threadNumber
        self.core_count: int = multiprocessing.cpu_count()
        
    def combine_fasta_files(self, fastaPath: str, decoyPath: str) -> None:
        self.logger.info(f'Combining target and decoy fasta files to {self.targetDecoyPath}')
        with open(fastaPath, 'r') as fasta, open(decoyPath, 'r') as decoy, open(self.targetDecoyPath, 'w') as output:
            for line in fasta:
                output.write(line)
            for line in decoy:
                output.write(line)
        
    def intergrate_filtered_psms_with_feature(self, baseName: str) -> pd.DataFrame:
        self.logger.info(f'Intergrating filtered psms with feature for {baseName}')
        target = pd.read_csv(f'{self.outputPath}/{baseName}/{baseName}_target_psms.tsv', sep='\t')
        decoy = pd.read_csv(f'{self.outputPath}/{baseName}/{baseName}_decoy_psms.tsv', sep='\t')
        pin = pd.read_csv(f'{self.outputPath}/{baseName}/{baseName}.pin', sep='\t')
        # for filtered PSMs
        filtered_target = target[target['q-value'] <= 0.01]
        filtered_decoy = decoy[decoy['q-value'] <= 0.01]
        psm = pd.concat([filtered_target, filtered_decoy], ignore_index=True)
        psm = psm[['PSMId', 'score', 'q-value', 'posterior_error_prob']]
        psm = pd.merge(psm, pin, left_on='PSMId', right_on='SpecId', how='left')
        psm['ScanNr'] = psm['ScanNr'].astype(int)
        psm.sort_values(by='ScanNr', inplace=True)
        psm.to_csv(f'{self.outputPath}/{baseName}/{baseName}_filtered_psms.tsv', sep='\t', index=False)
        # for philosopher input
        filtered_target = target[target['posterior_error_prob'] < 0.5]
        filtered_decoy = decoy[decoy['posterior_error_prob'] < 0.5]
        psm = pd.concat([filtered_target, filtered_decoy], ignore_index=True)
        psm = psm[['PSMId', 'score', 'q-value', 'posterior_error_prob']]
        psm = pd.merge(psm, pin, left_on='PSMId', right_on='SpecId', how='left')
        psm['ScanNr'] = psm['ScanNr'].astype(int)
        psm.sort_values(by='ScanNr', inplace=True)
        return psm
    
    def dataframe_to_pepxml(self, psm: pd.DataFrame, baseName: str) -> None:
        self.logger.info(f'Converting dataframe to pepxml for {baseName}')
        root = etree.Element("msms_pipeline_analysis")
        analysis_summary = etree.SubElement(root, "analysis_summary")
        msms_run_summary = etree.SubElement(root, "msms_run_summary")
        search_summary = etree.SubElement(msms_run_summary, "search_summary", precursor_mass_type="monoisotopic",
                                          fragment_mass_type="monoisotopic", search_engine="X! Tandem", search_engine_version="Sipros")
        search_database = etree.SubElement(search_summary, "search_database", local_path='targetDecoy.faa')
        for i, row in psm.iterrows():
            spectrum_query = etree.SubElement(msms_run_summary, "spectrum_query", 
                                            start_scan=str(row['ScanNr']), 
                                            assumed_charge=str(row['parentCharges']), 
                                            spectrum=str(row['SpecId']), 
                                            end_scan=str(row['ScanNr']), 
                                            index=str(i), 
                                            precursor_neutral_mass=str(row['ExpMass']), 
                                            retention_time_sec=str(row['retentiontime']))
            search_result = etree.SubElement(spectrum_query, "search_result")
            # split peptide sequence by "[" and "]" to get previous, current and next amino acids
            seqs = re.split(r'[\[\]]', row['Peptide']) 
            pros = row['Proteins'][1:-1].split(",")
            search_hit = etree.SubElement(search_result, "search_hit", 
                                        peptide=seqs[1], 
                                        massdiff=str(row['massErrors']), 
                                        calc_neutral_pep_mass=str(row['ExpMass']), 
                                        num_missed_cleavages=str(row['missCleavageSiteNumbers']), 
                                        num_tol_term="2", 
                                        protein_descr=pros[0], 
                                        num_tot_proteins=str(len(pros)), 
                                        hit_rank=str(row['ranks']), 
                                        protein=pros[0], 
                                        peptide_prev_aa=seqs[0], 
                                        peptide_next_aa=seqs[2],
                                        is_rejected="0")
            for k in range(1, len(pros)):
                alternative_protein = etree.SubElement(search_hit, "alternative_protein", 
                                                    protein_descr=pros[k], 
                                                    protein=pros[k], 
                                                    peptide_prev_aa=seqs[0], 
                                                    peptide_next_aa=seqs[2], 
                                                    num_tol_term="2")
            prob = str(1 - row['posterior_error_prob'])
            analysis_result = etree.SubElement(search_hit, "analysis_result", analysis="peptideprophet")
            # for 
            hyperscore = etree.SubElement(search_hit, "hyperscore", value=str(row['score']))
            nextscore = etree.SubElement(search_hit, "nextscore", value="0")
            expect = etree.SubElement(search_hit, "expect", value="0")
            peptideprophet_result = etree.SubElement(analysis_result, "peptideprophet_result", 
                                                    probability=prob, 
                                                    all_ntt_prob=f"({prob},{prob},{prob})")

        tree = etree.ElementTree(root)
        tree.write(f"{self.outputPath}/{baseName}/{baseName}.pep.xml", pretty_print=True, xml_declaration=True, encoding="UTF-8")
        
    def intergrate_and_convert(self, baseName: str) -> None:
        psm = self.intergrate_filtered_psms_with_feature(baseName)
        self.dataframe_to_pepxml(psm, baseName)

    def run_command(self, cmd:str, path:str = None) -> None:
        self.logger.info(f"Running command: {cmd}")
        original_cwd = os.getcwd()
        try:
            if path != None:
                os.chdir(path)
            output = subprocess.check_output(
                cmd, shell=True, stderr=subprocess.STDOUT)
            self.logger.info(output.decode())
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command execution failed: {e.output.decode()}")
            exit(1)
        finally:
            if path != None:
                os.chdir(original_cwd)  # Restore the original working directory

    def run(self) -> None:
        self.combine_fasta_files(self.fastaPath, self.decoyPath)
        threadNumber = self.core_count
        if self.threadNumber != None:
            threadNumber = self.threadNumber
        raw_file_parallel = min(10, threadNumber)
        
        self.logger.info(f'convert PSM tsv to pepXML with {raw_file_parallel} processes')
        # python function need processPoolExecutor to run in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            list(executor.map(self.intergrate_and_convert, self.baseNames)) # wait for all tasks to complete
        pepxmls_dir = os.path.join(self.outputPath, 'pepxmls')
        if not os.path.exists(pepxmls_dir):
            os.makedirs(pepxmls_dir)
        # Copy all generated .pepxml files to the pepxmls folder
        for baseName in self.baseNames:
            src = os.path.join(self.outputPath, baseName, f'{baseName}.pep.xml')
            dst = os.path.join(pepxmls_dir, f'{baseName}.pep.xml')
            self.logger.info(f'Copying {src} to {dst}')
            shutil.copy(src, dst)
        # infer proteins by all pepxml files
        pepxml_files = ' '.join([f'{baseName}/{baseName}.pep.xml' for baseName in self.baseNames])
        cmd = (
            f'{self.philosopherPath} workspace --init && '
            f'{self.philosopherPath} proteinprophet --maxppmdiff 2000000 --output combined {pepxml_files} && '
            f'{self.philosopherPath} database --annotate targetDecoy.faa --prefix Decoy_ && '
            f'{self.philosopherPath} filter --sequential --prot 0.01 --picked --tag Decoy_ --pepxml pepxmls --protxml combined.prot.xml --razor && '
            f'{self.philosopherPath} report'
        )
        self.run_command(cmd, self.outputPath)
        
        # filter and report for each raw file
        commands = []
        paths = []
        for baseName in self.baseNames:
            cmd = (
                f'{self.philosopherPath} workspace --init && '
                f'{self.philosopherPath} filter --sequential --prot 0.01 --picked --tag Decoy_ --dbbin .. --pepxml {baseName}.pep.xml --protxml ../combined.prot.xml --razor && '
                f'{self.philosopherPath} report'
            )
            commands.append(cmd)
            path=os.path.join(self.outputPath, baseName)
            paths.append(path)
        with concurrent.futures.ProcessPoolExecutor(max_workers=raw_file_parallel) as executor:
            list(executor.map(self.run_command, commands, paths))