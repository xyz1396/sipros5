from logging import Logger
import os
import shlex

from command_runner import run_logged_command
from thread_allocation import (
    available_cpu_count,
    effective_thread_count,
    thread_env_updates,
)


class filter:
    """Run the single Aerith post-search stage for the Python workflow.

    Aerith owns PSM filtering, FASTA chromatographic quantification, protein
    assembly, and protein-to-PSM reports. Python orchestrates the native
    commands; optional SIP negative-control filtering remains in the same
    native process and operates on Aerith's in-memory PSMs.
    """

    def __init__(self, baseNames: list[str], outputPath: str, aerithPath: str,
                 threadNumber: int, logger: Logger, decoyPrefix: str = "Decoy_",
                 ignorePCT: bool = False, dryrun: bool = False,
                 spectraPaths: list[str] | None = None,
                 fastaPath: str = "", decoyPath: str = "",
                 assembleProteins: bool = True, element: str = "",
                 negative_control: str = "", label_threshold: float = 2.0,
                 fixedCam: bool = True,
                 sipIsotope: str = "",
                 ptms: list[str] | None = None,
                 fixedPtms: list[str] | None = None,
                 maxPtmCount: int | None = None,
                 productTopIsotopes: int = 5,
                 quantTopIsotopes: int = 6) -> None:
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
        self.sipIsotope = sipIsotope
        self.ptms = [] if ptms is None else list(ptms)
        self.fixedPtms = [] if fixedPtms is None else list(fixedPtms)
        self.maxPtmCount = maxPtmCount
        self.productTopIsotopes = productTopIsotopes
        self.quantTopIsotopes = quantTopIsotopes
        if self.spectraPaths and len(self.spectraPaths) != len(self.baseNames):
            raise ValueError("Aerith requires one spectra path per sample")
        if self.assembleProteins and (not self.fastaPath or not self.decoyPath):
            raise ValueError(
                "Aerith protein assembly requires target and decoy FASTA paths"
            )
        if self.negative_control and not self.assembleProteins:
            raise ValueError(
                "Native negative-control filtering requires protein assembly"
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
        sip_isotope = getattr(self, "sipIsotope", "")
        if sip_isotope and sip_isotope != "R":
            arguments.extend(["--sip-isotope", sip_isotope])
        for fixed_ptm in getattr(self, "fixedPtms", []):
            arguments.extend(["--fixed-ptm", fixed_ptm])
        for ptm in getattr(self, "ptms", []):
            arguments.extend(["--ptm", ptm])
        max_ptm_count = getattr(self, "maxPtmCount", None)
        if max_ptm_count is not None:
            arguments.extend(["--max-ptm-count", str(max_ptm_count)])
        arguments.extend([
            "--product-top-isotopes",
            str(getattr(self, "productTopIsotopes", 5)),
            "--quant-top-isotopes",
            str(getattr(self, "quantTopIsotopes", 6)),
        ])
        negative_control = getattr(self, "negative_control", "")
        if negative_control:
            arguments.extend([
                "--negative-control", negative_control,
                "--label-threshold",
                str(getattr(self, "label_threshold", 2.0)),
            ])
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

    def postprocess_fasta_results(self) -> None:
        protein_report = os.path.join(
            self.outputPath, "combined_protein.tsv"
        )
        if not os.path.exists(protein_report):
            raise RuntimeError(
                f"Aerith did not create the native protein report: "
                f"{protein_report}"
            )
        if self.negative_control:
            for name in (
                "SIP_filtered_psms.tsv",
                "SIP_target_psms.tsv",
                "SIP_decoy_psms.tsv",
                "combined_protein_with_SIP_filtered_PSM.tsv",
            ):
                path = os.path.join(self.outputPath, name)
                if not os.path.exists(path):
                    raise RuntimeError(
                        "Aerith did not create native negative-control output: "
                        f"{path}"
                    )

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
