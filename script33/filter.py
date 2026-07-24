from logging import Logger
import os
import re
import shlex

from command_runner import run_logged_command
from thread_allocation import (
    available_cpu_count,
    configure_process_worker,
    effective_thread_count,
    thread_env_updates,
)

# This must run before pandas loads NumPy/OpenBLAS. Aerith receives the full
# native thread team; Python post-processing remains single-threaded.
configure_process_worker(1)

import pandas as pd


class filter:
    """Run the single Aerith post-search stage for the Python workflow.

    Aerith owns PSM filtering and, for regular FASTA searches, native
    chromatographic quantification and protein assembly. Python only prepares
    the command and optional denormalized convenience reports.
    """

    def __init__(self, baseNames: list[str], outputPath: str, aerithPath: str,
                 threadNumber: int, logger: Logger, decoyPrefix: str = "Decoy_",
                 ignorePCT: bool = False, dryrun: bool = False,
                 spectraPaths: list[str] | None = None,
                 fastaPath: str = "", decoyPath: str = "",
                 assembleProteins: bool = True, element: str = "",
                 negative_control: str = "", label_threshold: float = 2.0,
                 fixedCam: bool = True) -> None:
        self.aerithPath = aerithPath
        self.baseNames = baseNames
        self.outputPath = outputPath
        self.logger = logger
        self.core_count = available_cpu_count()
        self.threadNumber = effective_thread_count(threadNumber, self.core_count)
        self.decoyPrefix = decoyPrefix
        self.ignorePCT = ignorePCT
        self.dryrun = dryrun
        self.spectraPaths = [] if spectraPaths is None else list(spectraPaths)
        self.fastaPath = fastaPath
        self.decoyPath = decoyPath
        self.assembleProteins = assembleProteins
        self.element = element
        self.negative_control = negative_control
        self.label_threshold = label_threshold
        self.fixedCam = fixedCam
        if self.spectraPaths and len(self.spectraPaths) != len(self.baseNames):
            raise ValueError("Aerith requires one spectra path per sample")
        if self.assembleProteins and (not self.fastaPath or not self.decoyPath):
            raise ValueError(
                "Aerith protein assembly requires target and decoy FASTA paths"
            )

    def command(self) -> str:
        arguments = [self.aerithPath, "--decoy-prefix", self.decoyPrefix]
        if self.assembleProteins:
            arguments.extend([
                "--database", self.fastaPath,
                "--decoy-database", self.decoyPath,
                "--protein-output-dir", self.outputPath,
            ])
        else:
            arguments.extend(["--no-protein-assembly", "--filtered-only"])
        if self.ignorePCT:
            arguments.append("--ignore-pct")
        if not self.fixedCam:
            arguments.append("--no-fixed-cam")
        for index, baseName in enumerate(self.baseNames):
            sample = f"{self.outputPath}/{baseName}/{baseName}"
            arguments.extend([
                "--target-pin", f"{sample}_target.pin",
                "--decoy-pin", f"{sample}_decoy.pin",
                "--output-prefix", sample,
            ])
            if self.spectraPaths:
                arguments.extend(["--spectra", self.spectraPaths[index]])
        return shlex.join(arguments)

    def run_command(self, command: str) -> None:
        environment = thread_env_updates(self.threadNumber)
        environment["MKL_NUM_THREADS"] = str(self.threadNumber)
        run_logged_command(
            command,
            self.logger,
            env_updates=environment,
            cpu_cores=self.threadNumber,
        )

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

    def filter_sip_labeled_psms(
            self, filtered_psms: dict[str, pd.DataFrame]) -> pd.DataFrame:
        if not self.negative_control:
            self.logger.info("No negative control samples provided")
            return pd.DataFrame()
        negative_controls = [
            sample.strip() for sample in self.negative_control.split(",")
        ]
        missing = [
            sample for sample in negative_controls if sample not in self.baseNames
        ]
        if missing:
            self.logger.error(
                f"Negative control samples {missing} not found in baseNames "
                f"{self.baseNames}"
            )
            return pd.DataFrame()

        self.logger.info(
            f"Merging filtered SIP labeled PSMs with threshold "
            f"{self.label_threshold}"
        )
        all_psms = []
        for baseName in self.baseNames:
            sip_psms = filtered_psms[baseName].drop(
                columns=["score", "posterior_error_prob", "sqrtAbsDeltaRT"],
                errors="ignore",
            )
            sip_psms = sip_psms[
                (sip_psms["q-value"] <= 0.01)
                & (sip_psms["Label"] == 1)
                & (sip_psms["MS2IsotopicAbundances"] >= self.label_threshold)
            ].copy()
            sip_psms.drop(columns=["q-value"], inplace=True)
            if sip_psms.empty:
                self.logger.warning(
                    f"No SIP labeled PSMs passed first filters for {baseName}"
                )
                continue
            if baseName in negative_controls:
                sip_psms["Label"] = -1
                self.logger.info(
                    f"Added {len(sip_psms)} decoy SIP labeled PSMs from "
                    f"{baseName}"
                )
            else:
                self.logger.info(
                    f"Added {len(sip_psms)} target SIP labeled PSMs from "
                    f"{baseName}"
                )
            sip_psms["SampleName"] = baseName
            sip_psms["MS1IsotopicAbundances"] = sip_psms[
                "MS1IsotopicAbundances"
            ].clip(lower=0, upper=100)
            sip_psms.rename(columns={"PSMId": "SpecId"}, inplace=True)
            all_psms.append(sip_psms)
        if not all_psms:
            self.logger.warning(
                "No SIP labeled PSMs found that meet the filtering criteria"
            )
            return pd.DataFrame()

        merged = pd.concat(all_psms, ignore_index=True)
        input_path = f"{self.outputPath}/SIP.pin"
        merged.to_csv(input_path, sep="\t", index=False)
        self.logger.info(
            f"Merged {len(merged)} SIP labeled PSMs saved to {input_path}"
        )
        command = shlex.join([
            self.aerithPath,
            "--input", input_path,
            "--output-prefix", f"{self.outputPath}/SIP",
            "--no-protein-assembly",
            "--filtered-only",
            "--decoy-prefix", self.decoyPrefix,
        ])
        self.logger.info("Running Aerith filtering for SIP labeled PSMs")
        self.run_command(command)
        return pd.read_csv(
            f"{self.outputPath}/SIP_filtered_psms.tsv", sep="\t"
        )

    def match_psms_to_proteins(
            self, filtered_psms: dict[str, pd.DataFrame]) -> None:
        self.logger.info("Matching individual sample PSMs to proteins")
        protein_df = pd.read_csv(
            f"{self.outputPath}/combined_protein.tsv", sep="\t"
        )
        for baseName in self.baseNames:
            psm_df = filtered_psms[baseName].copy()
            psm_df = psm_df[psm_df["Label"] == 1]
            psm_df["PeptideSequence"] = psm_df["Peptide"].str.extract(
                r"\[([^\]]+)\]"
            )
            psm_df["Proteins"] = psm_df["Proteins"].str.strip("{}")
            psm_df = psm_df.assign(
                Protein=psm_df["Proteins"].str.split(",")
            ).explode("Protein")
            psm_df["Protein"] = (
                psm_df["Protein"].str.strip().map(self.sanitize_protein_id)
            )
            aggregate = psm_df.groupby("Protein").agg({
                "PeptideSequence": lambda values: ",".join(values.astype(str)),
                "MS1IsotopicAbundances": lambda values: ",".join(values.astype(str)),
                "MS2IsotopicAbundances": lambda values: ",".join(values.astype(str)),
                "log10_precursorIntensities": lambda values: ",".join(values.astype(str)),
            }).reset_index()
            prefix = f"{baseName}_"
            protein_df = protein_df.merge(
                aggregate.rename(columns={
                    "PeptideSequence": f"{prefix}PeptideSequences",
                    "MS1IsotopicAbundances": f"{prefix}MS1IsotopicAbundances",
                    "MS2IsotopicAbundances": f"{prefix}MS2IsotopicAbundances",
                    "log10_precursorIntensities":
                        f"{prefix}log10_precursorIntensities",
                }),
                on="Protein",
                how="left",
            )
            sample_proteins = pd.read_csv(
                f"{self.outputPath}/{baseName}/protein.tsv",
                sep="\t",
                usecols=["Protein", "Razor Intensity"],
            ).rename(columns={
                "Razor Intensity": f"{baseName}_ProteinAbundance",
            })
            protein_df = protein_df.merge(
                sample_proteins,
                on="Protein",
                how="left",
            )
            protein_df[f"{baseName}_ProteinAbundance"] = protein_df[
                f"{baseName}_ProteinAbundance"
            ].fillna(0.0)
        output_path = f"{self.outputPath}/combined_protein_with_PSM.tsv"
        protein_df.to_csv(output_path, sep="\t", index=False)
        self.logger.info(
            f"Updated protein information with filtered PSMs: {output_path}"
        )

    def match_sip_psms_to_proteins(self, sip_psms: pd.DataFrame) -> None:
        if sip_psms.empty:
            return
        self.logger.info("Matching SIP filtered PSMs to proteins")
        protein_df = pd.read_csv(
            f"{self.outputPath}/combined_protein.tsv", sep="\t"
        )
        sip_psms = sip_psms[sip_psms["Label"] == 1].copy()
        sip_psms["PeptideSequence"] = sip_psms["Peptide"].str.extract(
            r"\[([^\]]+)\]"
        )
        sip_psms["Proteins"] = sip_psms["Proteins"].str.strip("{}")
        sip_psms = sip_psms.assign(
            Protein=sip_psms["Proteins"].str.split(",")
        ).explode("Protein")
        sip_psms["Protein"] = (
            sip_psms["Protein"].str.strip().map(self.sanitize_protein_id)
        )
        for sample in sip_psms["SampleName"].unique():
            sample_psms = sip_psms[sip_psms["SampleName"] == sample]
            aggregate = sample_psms.groupby("Protein").agg({
                "PeptideSequence": lambda values: ",".join(values.astype(str)),
                "MS1IsotopicAbundances": lambda values: ",".join(values.astype(str)),
                "MS2IsotopicAbundances": lambda values: ",".join(values.astype(str)),
                "log10_precursorIntensities": lambda values: ",".join(values.astype(str)),
            }).reset_index().rename(columns={
                "PeptideSequence": f"SIP_{sample}_PeptideSequences",
                "MS1IsotopicAbundances": f"SIP_{sample}_MS1IsotopicAbundances",
                "MS2IsotopicAbundances": f"SIP_{sample}_MS2IsotopicAbundances",
                "log10_precursorIntensities":
                    f"SIP_{sample}_log10_precursorIntensities",
            })
            protein_df = protein_df.merge(aggregate, on="Protein", how="left")
        output_path = (
            f"{self.outputPath}/combined_protein_with_SIP_filtered_PSM.tsv"
        )
        protein_df.to_csv(output_path, sep="\t", index=False)
        self.logger.info(
            f"Updated protein information with SIP filtered PSMs: {output_path}"
        )

    def postprocess_fasta_results(self) -> None:
        filtered_psms = {
            baseName: pd.read_csv(
                f"{self.outputPath}/{baseName}/{baseName}_filtered_psms.tsv",
                sep="\t",
            )
            for baseName in self.baseNames
        }
        if all(psms.empty for psms in filtered_psms.values()):
            self.logger.warning(
                "No filtered PSMs found after Aerith; skipping protein annotation"
            )
            return
        protein_report = os.path.join(
            self.outputPath, "combined_protein.tsv"
        )
        if not os.path.exists(protein_report):
            raise RuntimeError(
                f"Aerith did not create the native protein report: "
                f"{protein_report}"
            )
        sip_psms = pd.DataFrame()
        if self.element and self.element != "R":
            sip_psms = self.filter_sip_labeled_psms(filtered_psms)
        self.match_psms_to_proteins(filtered_psms)
        if not sip_psms.empty:
            self.match_sip_psms_to_proteins(sip_psms)

    def run(self) -> None:
        if not self.baseNames:
            self.logger.info("Aerith: no jobs")
            return
        command = self.command()
        if self.assembleProteins and self.spectraPaths:
            operation = (
                "filtering, chromatographic quantification, and "
                "protein assembly"
            )
        elif self.assembleProteins:
            operation = "filtering and protein assembly"
        else:
            operation = "filtering only"
        self.logger.info(
            f"Running Aerith cross-sample {operation} with "
            f"{self.threadNumber} CPU cores"
        )
        if self.spectraPaths:
            self.logger.info(
                "Aerith DIA-NN device policy: CUDA preferred; "
                "automatic CPU fallback when CUDA is unavailable or fails; "
                "legacy Aerith RT modeling is skipped"
            )
        self.logger.info(command)
        if self.dryrun:
            return
        self.run_command(command)
        if self.assembleProteins:
            self.postprocess_fasta_results()
