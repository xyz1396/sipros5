from logging import Logger
import os
import concurrent.futures
import re
import shutil
import tempfile
from command_runner import run_logged_command
from thread_allocation import (
    ThreadAllocation,
    allocate_threads,
    available_cpu_count,
    configure_process_worker,
    effective_thread_count,
    thread_env_updates,
)

# This must run before pandas loads NumPy/OpenBLAS.  The Python process-pool
# phase is file-parallel, while external commands get explicit larger quotas.
configure_process_worker(1)

import pandas as pd
from lxml import etree


class assembly:
    def __init__(self, baseNames: list[str], philosopherPath:str, percolatorPath:str, fastaPath:str, decoyPath:str, outputPath: str,
                threadNumber: int, logger: Logger, element:str, negative_control = "", label_threshold = 2.0, decoyPrefix: str = "Decoy_") -> None:
        self.philosopherPath = philosopherPath
        self.percolatorPath = percolatorPath
        self.fastaPath = fastaPath
        self.decoyPath = decoyPath
        self.baseNames = baseNames
        self.outputPath = outputPath
        # merged fasta file for philosopher
        self.targetDecoyPath = f'{outputPath}/targetDecoy.faa'
        self.logger = logger
        self.negative_control = negative_control
        self.label_threshold = label_threshold
        self.core_count = available_cpu_count()
        self.threadNumber = effective_thread_count(threadNumber, self.core_count)
        self.element = element
        self.decoyPrefix = decoyPrefix
        
    def sanitize_protein_id(self, protein_id: str) -> str:
        protein_id = str(protein_id).strip()
        active_prefix = ""
        core = protein_id
        for prefix in (self.decoyPrefix, "DECOY_", "Decoy_"):
            if core.startswith(prefix):
                active_prefix = self.decoyPrefix
                core = core[len(prefix):]
                break
        if "|" not in core:
            safe = re.sub(r"\s+", "_", core)
            core = f"sp|{safe}|{safe}"
        return active_prefix + core

    def sanitize_fasta_header(self, line: str) -> str:
        header = line[1:].strip()
        if not header:
            return ">sp|UNKNOWN|UNKNOWN\n"
        parts = header.split(None, 1)
        protein_id = self.sanitize_protein_id(parts[0])
        suffix = f" {parts[1]}" if len(parts) > 1 else ""
        if " OS=" not in suffix:
            core = protein_id
            for prefix in (self.decoyPrefix, "DECOY_", "Decoy_"):
                if core.startswith(prefix):
                    core = core[len(prefix):]
                    break
            desc_parts = core.split("|")
            gene = desc_parts[2] if len(desc_parts) > 2 else desc_parts[-1]
            suffix += f" OS=Unknown OX=0 GN={gene} PE=4 SV=1"
        return f">{protein_id}{suffix}\n"

    def write_sanitized_fasta(self, source: str, output) -> bool:
        last_ends_with_newline = True
        with open(source, "r") as fasta:
            for line in fasta:
                out_line = self.sanitize_fasta_header(line) if line.startswith(">") else line
                output.write(out_line)
                last_ends_with_newline = out_line.endswith("\n")
        return last_ends_with_newline

    def combine_fasta_files(self, fastaPath: str, decoyPath: str) -> None:
        self.logger.info(f"Combining target and decoy fasta files to {self.targetDecoyPath}")
        with open(self.targetDecoyPath, "w") as output:
            target_ends_with_newline = self.write_sanitized_fasta(fastaPath, output)
            if not target_ends_with_newline:
                output.write("\n")
            self.write_sanitized_fasta(decoyPath, output)
        
    def intergrate_filtered_psms_with_feature(self, baseName: str) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        psm = psm.drop(columns=['SpecId'])
        psm['ScanNr'] = psm['ScanNr'].astype(int)
        psm.sort_values(by='ScanNr', inplace=True)
        psm.to_csv(f'{self.outputPath}/{baseName}/{baseName}_filtered_psms.tsv', sep='\t', index=False)
        # for philosopher input
        filtered_target = target[target['posterior_error_prob'] < 0.5]
        if self.decoyPrefix == 'DECOY_':
            filtered_decoy = decoy.iloc[0:0]
        else:
            filtered_decoy = decoy[decoy['posterior_error_prob'] < 0.5]
        psm_phi = pd.concat([filtered_target, filtered_decoy], ignore_index=True)
        psm_phi = psm_phi[['PSMId', 'score', 'q-value', 'posterior_error_prob']]
        psm_phi = pd.merge(psm_phi, pin, left_on='PSMId', right_on='SpecId', how='left')
        psm_phi['ScanNr'] = psm_phi['ScanNr'].astype(int)
        psm_phi.sort_values(by='ScanNr', inplace=True)
        return psm, psm_phi
    
    def dataframe_to_pepxml(self, psm: pd.DataFrame, baseName: str) -> None:
        self.logger.info(f'Converting dataframe to pepxml for {baseName}')
        root = etree.Element("msms_pipeline_analysis")
        analysis_summary = etree.SubElement(root, "analysis_summary")
        msms_run_summary = etree.SubElement(root, "msms_run_summary")
        search_summary = etree.SubElement(msms_run_summary, "search_summary", precursor_mass_type="monoisotopic",
                                          fragment_mass_type="monoisotopic", search_engine="X! Tandem", search_engine_version="Sipros")
        outputPath = os.path.abspath(self.outputPath)
        search_database = etree.SubElement(search_summary, "search_database", local_path=f'{outputPath}/targetDecoy.faa', type="AA")
        parameter = etree.SubElement(search_summary, "parameter", name="database_name", value=f'{outputPath}/targetDecoy.faa')
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
            if len(seqs) < 3:
                raise ValueError(f"Unexpected peptide format: {row['Peptide']}")
            # handle terminal amino acids to avoid empty strings making crash in philosopher
            if not seqs[0] or seqs[0].strip() == "":
                seqs[0] = "-"
            if not seqs[2] or seqs[2].strip() == "":
                seqs[2] = "-"
            seq = seqs[1]
            # remove PTM in the pep seq
            seq = re.sub(r'[^a-zA-Z]', '', seq)
            pros = [self.sanitize_protein_id(pro.strip()) for pro in row['Proteins'][1:-1].split(",") if pro.strip()]
            if not pros:
                continue
            search_hit = etree.SubElement(search_result, "search_hit", 
                                        peptide=seq, 
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
            hyperscore = etree.SubElement(search_hit, "hyperscore", value=str(row['score']))
            nextscore = etree.SubElement(search_hit, "nextscore", value="0")
            expect = etree.SubElement(search_hit, "expect", value="0")
            peptideprophet_result = etree.SubElement(analysis_result, "peptideprophet_result", 
                                                    probability=prob, 
                                                    all_ntt_prob=f"({prob},{prob},{prob})")

        tree = etree.ElementTree(root)
        tree.write(f"{self.outputPath}/{baseName}/{baseName}.pep.xml", pretty_print=True, xml_declaration=True, encoding="UTF-8")
        
    def intergrate_and_convert(self, baseName: str) -> tuple[str, pd.DataFrame]:
        psm, psm_phi = self.intergrate_filtered_psms_with_feature(baseName)
        self.dataframe_to_pepxml(psm_phi, baseName)
        return baseName, psm
        
    def filterSIPlabeledPSMs(self, filteredPSMsDict: dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Merge all filtered PSMs with SIP labeling information.
        PSMs with MS2abundance > label_threshold are kept.
        PSMs from negative control samples are marked as decoy.
        """
        # Parse negative control sample names if provided
        negative_controls = []
        if self.negative_control != "" and self.negative_control != None:
            negative_controls = [sample.strip() for sample in self.negative_control.split(',')]
            if not all(control in self.baseNames for control in negative_controls):
                missing_controls = [control for control in negative_controls if control not in self.baseNames]
                self.logger.error(f'Negative control samples {missing_controls} not found in baseNames {self.baseNames}')
                return pd.DataFrame()
        else:
            self.logger.info('No negative control samples provided')
            return pd.DataFrame()

        self.logger.info(f'Merging filtered SIP labeled PSMs with threshold {self.label_threshold}')
        all_psms = []
        
        # Process each sample
        for baseName in self.baseNames:
            SIP_df = filteredPSMsDict[baseName].drop(columns=['score','posterior_error_prob'])
            SIP_df = SIP_df[(SIP_df['q-value'] <= 0.01) & (SIP_df['Label'] == 1) &
                            (SIP_df['MS2IsotopicAbundances'] >= self.label_threshold)]
            SIP_df.drop(columns=['q-value'], inplace=True)
            if SIP_df.empty:
                self.logger.warning(f'No SIP labeled PSMs passed first filters for {baseName}')
                continue
            # Mark PSMs from negative controls as decoy
            if baseName in negative_controls:
                SIP_df['Label'] = -1
                self.logger.info(f'Added {len(SIP_df)} decoy SIP labeled PSMs from {baseName}')
            else:
                self.logger.info(f'Added {len(SIP_df)} target SIP labeled PSMs from {baseName}')
            SIP_df.rename(columns={'PSMId': 'SpecId'}, inplace=True)  
            all_psms.append(SIP_df)
        merged_psms = pd.DataFrame()
        # Merge all PSMs
        if all_psms:
            merged_psms = pd.concat(all_psms, ignore_index=True)
            output_path = f'{self.outputPath}/SIP.pin'
            merged_psms.to_csv(output_path, sep='\t', index=False)
            self.logger.info(f'Merged {len(merged_psms)} SIP labeled PSMs saved to {output_path}')
        else:
            self.logger.warning('No SIP labeled PSMs found that meet the filtering criteria')
            return pd.DataFrame()
        cmd = f'{self.percolatorPath} --only-psms --no-terminate -Y --num-threads {self.threadNumber} \
                    --results-psms {self.outputPath}/SIP_target_psms.tsv \
                    --decoy-results-psms {self.outputPath}/SIP_decoy_psms.tsv \
                    {self.outputPath}/SIP.pin'
        self.logger.info(f'Running percolator for SIP labeled PSMs')
        self.run_command(cmd)
        self.logger.info(f'Intergrating filtered SIP labeled PSMs with feature')
        target = pd.read_csv(f'{self.outputPath}/SIP_target_psms.tsv', sep='\t')
        decoy = pd.read_csv(f'{self.outputPath}/SIP_decoy_psms.tsv', sep='\t')
        # for filtered SIP labeled PSMs
        filtered_target = target[target['q-value'] <= 0.01]
        filtered_decoy = decoy[decoy['q-value'] <= 0.01]
        psm = pd.concat([filtered_target, filtered_decoy], ignore_index=True)
        psm = psm[['PSMId', 'score', 'q-value', 'posterior_error_prob']]
        psm = pd.merge(psm, merged_psms, left_on='PSMId', right_on='SpecId', how='left')
        # Split baseName from PSMId column
        psm['SampleName'] = psm['PSMId'].str.split('.').str[0]
        psm.drop(columns=['SpecId'], inplace=True)
        psm['ScanNr'] = psm['ScanNr'].astype(int)
        psm.sort_values(by=['SampleName', 'ScanNr'], inplace=True)
        # In case Non-carbon elements SIP exceed 100. Clamp to [0, 100]
        psm['MS1IsotopicAbundances'] = psm['MS1IsotopicAbundances'].clip(lower=0, upper=100)
        psm.to_csv(f'{self.outputPath}/SIP_filtered_psms.tsv', sep='\t', index=False)
        return psm

    def match_PSMs_to_proteins(self, filteredPSMsDict: dict[str, pd.DataFrame]) -> None:
        """
        Add peptide information from individual sample with PSMs passed decoy filtering to protein.tsv.
        For each sample and each protein, add columns for peptideSequence, MS1IsotopicAbundances, 
        MS2IsotopicAbundances, and log10_precursorIntensities from matched PSMs.
        Multiple values are separated by commas.
        """
        self.logger.info('Matching individual sample with PSMs to proteins')
        protein_df = pd.read_csv(f'{self.outputPath}/protein.tsv', sep='\t')
        for baseName in self.baseNames:
            if baseName in filteredPSMsDict:
                psm_df = filteredPSMsDict[baseName].copy()
                psm_df = psm_df[psm_df['Label'] == 1]
                psm_df['PeptideSequence'] = psm_df['Peptide'].str.extract(r'\[([^\]]+)\]')
                # Explode proteins
                psm_df['Proteins'] = psm_df['Proteins'].str.strip('{}')
                psm_df = psm_df.assign(Protein=psm_df['Proteins'].str.split(',')).explode('Protein')
                psm_df['Protein'] = psm_df['Protein'].str.strip().map(self.sanitize_protein_id)
                # Group by protein
                agg_df = psm_df.groupby('Protein').agg({
                    'PeptideSequence': lambda x: ','.join(x.astype(str)),
                    'MS1IsotopicAbundances': lambda x: ','.join(x.astype(str)),
                    'MS2IsotopicAbundances': lambda x: ','.join(x.astype(str)),
                    'log10_precursorIntensities': lambda x: ','.join(x.astype(str)),
                }).reset_index()
                col_prefix = f'{baseName}_'
                # Merge with protein_df
                protein_df = protein_df.merge(
                    agg_df.rename(columns={
                        'PeptideSequence': f'{col_prefix}PeptideSequences',
                        'MS1IsotopicAbundances': f'{col_prefix}MS1IsotopicAbundances',
                        'MS2IsotopicAbundances': f'{col_prefix}MS2IsotopicAbundances',
                        'log10_precursorIntensities': f'{col_prefix}log10_precursorIntensities',
                    }),
                    on='Protein', how='left'
                )
            else:
                self.logger.warning(f'No filtered PSMs found for {baseName}, skipping protein matching')
        output_path = f'{self.outputPath}/protein_with_PSM.tsv'
        protein_df.to_csv(output_path, sep='\t', index=False)
        self.logger.info(f'Updated protein information with PSMs passed decoy filtering saved to {output_path}')

    def match_SIP_filtered_PSMs_to_proteins(self, SIPfilteredPSM: pd.DataFrame) -> None:
        """
        Add peptide information with SIP filtered PSMs to protein.tsv.
        For each SIP sample and each protein, add columns for peptideSequence, MS1IsotopicAbundances, 
        MS2IsotopicAbundances, and log10_precursorIntensities from matched SIP PSMs.
        Multiple values are separated by commas.
        """
        self.logger.info('Matching individual sample with SIP filtered PSMs to proteins')
        protein_df = pd.read_csv(f'{self.outputPath}/protein.tsv', sep='\t')  
        if SIPfilteredPSM is not None and not SIPfilteredPSM.empty:
            SIPfilteredPSM = SIPfilteredPSM[SIPfilteredPSM['Label'] == 1].copy()
            SIPfilteredPSM['PeptideSequence'] = SIPfilteredPSM['Peptide'].str.extract(r'\[([^\]]+)\]')
            SIPfilteredPSM['Proteins'] = SIPfilteredPSM['Proteins'].str.strip('{}')
            SIPfilteredPSM = SIPfilteredPSM.assign(Protein=SIPfilteredPSM['Proteins'].str.split(',')).explode('Protein')
            SIPfilteredPSM['Protein'] = SIPfilteredPSM['Protein'].str.strip().map(self.sanitize_protein_id)
            sip_samples = SIPfilteredPSM['SampleName'].unique()
            # Check if SIP samples are in baseNames
            if not all(sample in self.baseNames for sample in sip_samples):
                missing_samples = [sample for sample in sip_samples if sample not in self.baseNames]
                self.logger.warning(f'SIP samples {missing_samples} not found in baseNames {self.baseNames}')            
            for sample in sip_samples:
                sample_psms: pd.DataFrame = SIPfilteredPSM[SIPfilteredPSM['SampleName'] == sample]
                agg_df = sample_psms.groupby('Protein').agg({
                    'PeptideSequence': lambda x: ','.join(x.astype(str)),
                    'MS1IsotopicAbundances': lambda x: ','.join(x.astype(str)),
                    'MS2IsotopicAbundances': lambda x: ','.join(x.astype(str)),
                    'log10_precursorIntensities': lambda x: ','.join(x.astype(str)),
                }).reset_index()
                sip_col_prefix = f'SIP_{sample}_'
                agg_df = agg_df.rename(columns={
                    'PeptideSequence': f'{sip_col_prefix}PeptideSequences',
                    'MS1IsotopicAbundances': f'{sip_col_prefix}MS1IsotopicAbundances',
                    'MS2IsotopicAbundances': f'{sip_col_prefix}MS2IsotopicAbundances',
                    'log10_precursorIntensities': f'{sip_col_prefix}log10_precursorIntensities',
                })
                protein_df = protein_df.merge(agg_df, on='Protein', how='left')
        else:
            self.logger.warning('No SIP filtered PSMs found, skipping SIP protein matching')
            return
        # Save the updated protein dataframe
        output_path = f'{self.outputPath}/protein_with_SIP_filtered_PSM.tsv'
        protein_df.to_csv(output_path, sep='\t', index=False)
        self.logger.info(f'Updated protein information with SIP filtered PSM saved to {output_path}')    

    def run_command(self, cmd: str, path: str = "",
                    threads: int | None = None) -> None:
        thread_count = self.threadNumber if threads is None else threads
        run_logged_command(
            cmd,
            self.logger,
            cwd=path or None,
            env_updates=thread_env_updates(thread_count),
            cpu_cores=thread_count,
        )

    def run_parallel_commands(self, commands: list[str], paths: list[str],
                              phase: str) -> ThreadAllocation:
        if len(commands) != len(paths):
            raise ValueError('commands and paths must have the same length')
        allocation = allocate_threads(self.threadNumber, len(commands))
        if allocation.worker_count == 0:
            self.logger.info(f'{phase}: no jobs')
            return allocation
        minimum = min(allocation.task_threads)
        maximum = max(allocation.task_threads)
        per_job = str(minimum) if minimum == maximum else f'{minimum}-{maximum}'
        unit = 'core' if minimum == maximum == 1 else 'cores'
        self.logger.info(
            f'{phase}: up to {allocation.worker_count} concurrent processes; '
            f'{per_job} CPU {unit} per process; '
            f'{allocation.peak_threads}/{self.threadNumber} cores allocated at peak'
        )
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=allocation.worker_count) as executor:
            list(executor.map(
                self.run_command, commands, paths, allocation.task_threads
            ))
        return allocation

    def run(self) -> None:
        self.combine_fasta_files(self.fastaPath, self.decoyPath)

        conversion_allocation = allocate_threads(
            self.threadNumber, len(self.baseNames)
        )
        filteredPSMsDict = {}
        if conversion_allocation.worker_count == 0:
            self.logger.warning('No samples available for protein assembly')
            return
        self.logger.info(
            f'Convert PSM TSV to pepXML with up to '
            f'{conversion_allocation.worker_count} single-threaded processes'
        )
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=conversion_allocation.worker_count,
                initializer=configure_process_worker,
                initargs=(1,)) as executor:
            future_to_baseName = {executor.submit(self.intergrate_and_convert, baseName): baseName for baseName in self.baseNames}
            for future in concurrent.futures.as_completed(future_to_baseName):
                baseName, psm_df = future.result()
                filteredPSMsDict[baseName] = psm_df
        # filter SIP labeled PSMs for SIP search
        SIPfilteredPSMs = pd.DataFrame()
        if (self.element != None and self.element != ""):
            SIPfilteredPSMs = self.filterSIPlabeledPSMs(filteredPSMsDict)
        if all(psm_df.empty for psm_df in filteredPSMsDict.values()):
            self.logger.warning("No filtered PSMs found after Percolator; skipping protein inference")
            return
        
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
        # TEMP and TMP are for windows
        tmpdir = os.environ.get('TMPDIR') or os.environ.get('TEMP') or os.environ.get('TMP') or tempfile.gettempdir()
        tmpdir = os.path.join(tmpdir, 'philosopher_tmp')
        self.logger.info(f'Inferring proteins by all pepxml files in {pepxmls_dir} with tmpdir {tmpdir}')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        cmd = (
            f'{self.philosopherPath} workspace --clean --nocheck && '
            f'{self.philosopherPath} workspace --init --nocheck --temp {tmpdir} && '
            f'{self.philosopherPath} proteinprophet --maxppmdiff 2000000 --output combined {pepxml_files} && '
            f'{self.philosopherPath} database --annotate targetDecoy.faa --prefix {self.decoyPrefix} && '
            f'{self.philosopherPath} filter --sequential --prot 0.01 --picked --tag {self.decoyPrefix} --pepxml pepxmls --protxml combined.prot.xml --razor && '
            f'{self.philosopherPath} report'
        )
        self.run_command(cmd, self.outputPath, self.threadNumber)
        
        # match decoy filtered PSMs to proteins
        self.match_PSMs_to_proteins(filteredPSMsDict)
        # match SIP filtered PSMs to proteins for SIP search
        if (self.element != None and self.element != "" and self.negative_control != "" and not SIPfilteredPSMs.empty):
            self.match_SIP_filtered_PSMs_to_proteins(SIPfilteredPSMs)
        
        # filter and report for each raw file
        commands = []
        paths = []
        # in case tmpdir was deleted
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
            self.logger.info(f'Re-created tmpdir {tmpdir}')
        for baseName in self.baseNames:
            tmp_subdir = tempfile.mkdtemp(prefix=f"{baseName}_", dir=tmpdir)
            cmd = (
                f'{self.philosopherPath} workspace --clean --nocheck && '
                f'{self.philosopherPath} workspace --init --nocheck --temp {tmp_subdir} && '
                f'{self.philosopherPath} filter --sequential --prot 0.01 --picked --tag {self.decoyPrefix} --dbbin .. --pepxml {baseName}.pep.xml --protxml ../combined.prot.xml --razor && '
                f'{self.philosopherPath} report'
            )
            commands.append(cmd)
            path=os.path.join(self.outputPath, baseName)
            paths.append(path)
        self.run_parallel_commands(
            commands, paths, 'Per-sample Philosopher report'
        )
        # remove tmpdir
        shutil.rmtree(tmpdir, ignore_errors=True)
