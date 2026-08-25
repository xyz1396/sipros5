<p align="center">
  <img src="wf33/sipros_logo.png" alt="Sipros5 logo" width="240">
</p>

## Sipros5 Setup Guide

### 1. Create Conda Environment

```bash
conda install bioconda::sipros
```

### 2. Download Raw Files

```bash
mkdir raw
# Download raw file with 1% 13C
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD041414/Pan_062822_X1iso5.raw -P raw
# Download raw file with 50% 13C  
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD041414/Pan_052322_X13.raw -P raw
```

### 3. Download E. coli Protein FASTA Sequence

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz
gunzip UP000000625_83333.fasta.gz -c > Ecoli.fasta
```

### 4. Example Commands

#### Regular Search

```bash
siproswf --regular-fasta-search -i raw/Pan_062822_X1iso5.raw -f Ecoli.fasta -o regular_output
```

#### Extract protein sequences identified in Regular search

- This step is particularly useful when your protein FASTA is large (for example, several GB in metaproteomics studies). The `regular_output/combined_protein.tsv` file can be replaced with a protein report from another search engine as long as the first column contains the protein identifier.
- If you are working with a small FASTA, you can skip this extraction step and use the original FASTA for the label search.

```bash
extractPro Ecoli.fasta regular_output/combined_protein.tsv db.faa
```

#### Label Search

```bash
siproswf --sip-fasta-search -i raw -f db.faa -e C13 -o sip_output
```

#### Label Search with negative control using unlabeled sample

```bash
siproswf --sip-fasta-search -i raw -f db.faa -e C13 --negative_control Pan_062822_X1iso5 -o sip2_output
```

Use a comma-separated `-r` list to search disjoint SIP enrichments in one
joint filtering run, for example `-r 1,2,3,49,50,51 -p 1`.

Run `siproswf` without arguments (or use `siproswf --gui`) to open the native
Dear ImGui interface. The GUI provides Regular FASTA, SIP FASTA, and Fast SIP
modes, advanced search settings, command preview, live logs, and cancellation.
Supplying workflow arguments keeps execution headless for scripts and batch
systems. Existing commands that infer regular versus SIP search from `-e`
remain supported.

#### Search chemistry and defaults

The Regular and SIP search profiles, residue/PTM chemistry, isotope
distributions, digestion rules, and default tolerances live in
`sipros/src/proNovoConfig.cpp`. Runtime inputs such as the FASTA database,
mass-tolerance overrides, SIP isotope, abundance range, and abundance step are
passed explicitly on the command line. Missing SIP controls are rejected.

Regular FASTA search enables oxidation (`~`, M) by default, matching the
validated benchmark's residue-level variable-modification space, with at most
three variable PTMs per peptide. Deamidation (`!`, N/Q) remains
available explicitly. Use a
repeatable `--ptm` option to replace that default set with compiled PTMs, and
`--max-ptm-count` to change the per-peptide limit. The `default` selector adds
the Regular defaults to an explicit list; `none` disables all variable PTMs;
and `all` selects every compatible variable PTM. Descriptive names avoid shell
quoting problems with symbols such as `>`, `&`, and `$`.

Regular HDF5 search derives candidates only from the acquisition isolation
window. For charges 1-4, every theoretical precursor inside that window enters
a fast fragment-index pass; MS1 peaks never create or prioritize candidates.
The gate uses the validated yeast DDA settings: the 200 most intense processed
peaks, singly charged b/y ions, and
at least four fragment matches. Exact MVH, XCorr, and WDP scoring follows for
gate survivors. The usual single acquisition window stays as a peptide-mass
range until the gate passes, avoiding a large temporary candidate vector.

After scoring and final per-scan pruning, each retained PSM is independently
matched against the linked parent MS1 acquisition plus or minus five MS1
acquisitions. The closest match by absolute retention-time difference is used,
with mass error and intensity as tie-breakers. PIN output records that distance
in seconds as `absPrecursorRtDiffSeconds`, immediately after
`isotopicPeakNumbers`; `-1` means no peptide-specific precursor peak was
matched. In that case, the linked parent MS1 acquisition and the isolation
window center anchor MS1-abundance imputation while the distance remains `-1`.
The imputed abundance is retained, but `isotopicPeakNumbers` and
`MS1IsotopeFitScore` are forced to zero because no peptide-specific precursor
was matched. If the center cannot anchor an observed MS1 envelope, the
theoretical isotope m/z closest to the isolation-window center supplies a
nominal-shift abundance estimate corrected for the peptide's natural-isotope
background. `isotopicAbundanceDiffs` remains the signed
`MS1IsotopicAbundances - MS2IsotopicAbundances` value. Up to five greedy
residual-spectrum winners are also retained through `ddaResidualRank` and
`ddaResidualScore`; each winner's matched fragment peaks are removed before
selecting the next one. Isolation-window discovery is the only Regular
candidate path; there is no precursor-source compatibility switch.

Carbamidomethyl-Cys is the default fixed PTM. A repeatable `--fixed-ptm`
option replaces the fixed default set: use `--fixed-ptm none` for natural Cys,
or `--fixed-ptm carbamidomethyl` to open it explicitly.

```bash
# Default: oxidation
siproswf -i sample.h5 -f proteins.fasta -o regular_output

# Oxidation plus explicit deamidation
siproswf -i sample.h5 -f proteins.fasta -o deamidation_output \
  --ptm default --ptm deamidation

# Defaults plus phosphorylation and acetylation
siproswf -i sample.h5 -f proteins.fasta -o ptm_output \
  --ptm default --ptm phosphorylation --ptm acetylation \
  --max-ptm-count 3

# Natural Cys: close the default fixed carbamidomethyl PTM
siproswf -i sample.h5 -f proteins.fasta -o natural_cys_output \
  --fixed-ptm none

# Show every compiled selector, token, and site
bin/sipros search-fasta --list-ptms
```

The restored catalog contains phosphorylation and its two deterministic
neutral-loss forms, acetylation, mono/di/trimethylation, S-nitrosylation,
nitration, and beta-methylthiolation in addition to the two defaults. IAA
carbamidomethylation uses the `/` token when natural Cys is selected; it is
excluded from variable search while the equivalent fixed PTM is open.
Bracketed absolute or delta modification masses are translated through this
same catalog before theoretical or experimental spectrum masses are
calculated.

Elemental formulas also retain isotope-source provenance. Amino-acid atoms are
`Biosynthetic` and follow the selected SIP abundance; carbamidomethyl atoms are
`ReagentNatural`, and peptide-terminal H2O is `DigestionSolvent`. Both natural
sources always use the natural isotope distribution of the real CHONPS element.
Generated spectra libraries record this contract in `chemistry_profile_id`, and
`search-spectra` rejects libraries made with a different or legacy chemistry.
Experimental spectra generation defaults to SIP abundances `0-100` at one
percentage-point intervals. It writes one memory-mapped `spectra.sfi` target
index and, with `--decoy`, one `spectra_Decoy.sfi` generated-decoy index; each
record stores its own abundance and retention time. Input is strictly Aerith
`*_filtered_psms.tsv`; legacy `psm.tsv` and missing confidence-column fallbacks
are not supported. For every modified-peptide mass class and precursor charge,
the representative PSM is selected by lowest `posterior_error_prob`, then
highest `SVMscore`, then highest `WDPscores`. SFI v6 stores only the top
three peaks by default in every precursor and product-ion isotope envelope
(`--envelope-top-n` / workflow `--sfi-envelope-top-n`), and stores fragment m/z
at 0.001-Da fixed-point resolution in an 8-byte hot fragment record. Fragment
position and b/y kind share otherwise unused bits, while the less frequently
read experimental intensities use a sparse sidecar. The index uses a sorted
4-byte packed product-ion index partitioned into sparse five-minute RT segments.
Candidate gating reads only segments overlapping the configured RT window.
Spectrum generation, precursor ordering, compact-array flattening,
product-index blocks, payload validation, and SFI publication are parallel;
generation logs report the worker count, fragment/posting counts, index-build
time, checksum/layout time, parallel-write time, and final GiB. HDF5
theoretical-spectrum libraries are not supported. Spectra search runs the target and decoy indexes
in parallel. It creates charge-1--4 candidates only from the acquisition
isolation window plus the SFI RT range, then applies the same top-200,
singly-charged b/y, four-product-ion first gate as optimized regular search.
No MS1 peak is required before MVH, Xcorr, and WDP scoring. For retained PSMs,
it searches the parent MS1 +/-5 acquisitions and exposes the absolute
retention-time distance in seconds (or -1 when unmatched) as a model feature.
It collects gate survivors in per-thread record bitsets and performs the final
unique-record reduction in parallel, avoiding a serial candidate-ID sort; the
per-file timing table reports measured wall time, process CPU time, and
parallel speedup for each search and output stage, followed by exact gate,
cascade, precursor-match, and PIN-row counts. It precomputes each surviving
SFI record's compact MVH ion list once, retains the
best 150 MVH candidates per scan by default for the Xcorr/WDP cascade
(`--mvh-cascade-top-n` changes this limit), and
materializes compact spectra only for that bounded cascade. Final WDP scoring
regenerates the full high-resolution b/y isotope envelopes from the peptide,
per-record SIP abundance, and SFI chemistry metadata just as regular search
does; it does not rank from the truncated SFI envelope. Matched-ion and
spectrum-shape features consume those same full WDP envelopes. The precursor
feature envelope is reconstructed with the established `b_(n-1) * y1` convolution
and retained after the much larger product envelopes are released. It writes one
`_target.pin` and one `_decoy.pin` per sample.

The workflow's fast SIP mode performs the complete two-stage search in one
command. It first runs the regular target/decoy FASTA cache search and Aerith
PSM filtering for every sample. It then uses only accepted target rows from
the resulting `*_filtered_psms.tsv` files to build one target and one generated
decoy SFI covering the requested SIP range, searches both indexes, and runs a
final cross-sample Aerith report. Existing HDF5 inputs are linked into the
regular stage instead of copied, and the same scan files are reused by spectra
generation, spectra search, and final feature calculation.
The workflow divides a CPU-thread budget among concurrent target/decoy jobs;
the configured thread limits are passed to each worker runtime.

```bash
siproswf -i 'T01.h5,T02.h5,T03.h5,X1.h5,X2.h5,X3.h5' \
  -f yeast.fasta --fast-sip-search -e C13 -r 0-100 -p 1 -t 48 \
  --aerith-sample-parallelism 3 \
  --negative_control X1,X2,X3 -o fast_sip_output
```

Regular-search reports, including peptide and protein TSVs, are written below
`regular/`; SIP spectra-search reports are written below `spectra_search/`.
The single target and generated-decoy SFI files are also written directly in
`spectra_search/`. Fast mode keeps a keyed cache of the regular-search DIA-NN
spectrum and RT predictions in `regular/`; regular-search decoy predictions
are not persisted. Aerith writes every unique target PIN peptide-charge form
to that cache during the regular-filtering prediction pass, so cache population
does not invoke the spectrum or RT model a second time. Regular cache
preparation writes `target.sfi` first. When
the sibling `decoy.sfi` is built, Sipros automatically loads `target.sfi` and
omits every canonical naked-peptide overlap (with I/J/L equivalence) before
building precursor and fragment postings. The target-index identity is part of
the decoy cache fingerprint, so a stale unguarded decoy cache is rebuilt. The
final Aerith pass reuses target predictions and predicts generated spectra-search
decoys in memory without persisting them or writing prediction catalog TSVs.
`--aerith-sample-parallelism` controls how many
samples Aerith processes concurrently; larger values can use more RAM. The
workflow logs the accepted regular-search target PSM count for each sample and
emits a warning when a sample has zero accepted target PSMs.

FASTA precursor estimates sum the expected exact isotope-mass shifts from all
source-aware CHONPS atoms, divide by the expected nominal-neutron shift to get
a peptide-specific spacing, and round the combined nominal shift once. The
same composition-weighted spacing defines that peptide's configured precursor
isotope windows; FASTA assignment does not fall back to a global target-isotope
spacing.

SIP FASTA search evaluates only precursor hypotheses stored by Raxport in each
HDF5 scan. It does not synthesize a precursor candidate from the reaction
precursor when the candidate list is empty.

### 5. Output Files

- `SIP_filtered_psms.tsv`: target PSMs from all non-control samples that pass the unlabeled negative-control filter (1% FDR). It reports the final `SVMscore`, the compact search `Peptide`, a human-readable `ModifiedPeptide`, `AssignedModifications`, and SIP element labeling percentages (`MS1IsotopicAbundances`, `MS2IsotopicAbundances`); training-only `Label` and `diffScores` fields are not exported. `isotopicPeakNumbers` is the raw number of extracted MS1 isotope peaks for a peptide-specific matched precursor and is zero for isolation-center imputation. `MS1IsotopeFitScore` is the theoretical-envelope coverage from `0` to `1`, set to `0` unless at least two compatible peaks and a peptide-specific precursor match are present; scores of at least `0.02` pass the MS1 abundance-fit validity threshold. MS1IsotopicAbundances are more sensitive; MS2IsotopicAbundances are more accurate.
- `SIP_target_psms.tsv`, `SIP_decoy_psms.tsv`: all target and negative-control candidates used by the secondary SIP SVM, written in the same `PSMId`, `SVMscore`, `q-value`, `posterior_error_prob`, `peptide`, `modifiedPeptide`, `assignedModifications`, and `proteinIds` format as the corresponding sample-subdirectory score tables.
- `combined_protein_with_SIP_filtered_PSM.tsv`: maps unlabeled negative-control filtered PSMs to the proteins identified in each sample.
- For each raw-file subdirectory:
  - `<sample>.h5`: Raxport scan data.
  - `<sample>_target.pin`, `<sample>_decoy.pin`: target and decoy search intermediates.
  - `<sample>_target_psms.tsv`, `<sample>_decoy_psms.tsv`: Aerith reranked `SVMscore` tables with compact and human-readable modified peptide forms and explicit assigned modifications.
  - `<sample>_filtered_psms.tsv`: Aerith 1% FDR PSMs with original search and RT features, `ModifiedPeptide`, and `AssignedModifications`.
  - `psm.tsv`: accepted PSM report with precursor
    intensity, assigned modifications, search class, protein coordinates, and
    FASTA annotations.
  - `ion.tsv`, `modified_peptide.tsv`, `peptide.tsv`: accepted PSMs aggregated
    by peptide ion, modified peptide, and naked peptide, with probabilities,
    spectral counts, precursor intensities, modifications, and protein mappings.
  - `protein.tsv`, `protein.fas`: Aerith's per-sample picked-FDR/razor protein report and identified FASTA entries.
  - `*_filtered_psms.tsv`: PSMs passing 1% FDR decoy filtering with `isotopicPeakNumbers`, `MS1IsotopeFitScore`, `MS1IsotopicAbundances`, `MS2IsotopicAbundances`, and, for regular and SIP FASTA search, `unweighted_spectral_entropy`, `delta_RT_loess_real`, and `pred_RT_real_units`.
- At the workflow root, `combined_psm.tsv`, `combined_ion.tsv`,
  `combined_modified_peptide.tsv`, `combined_peptide.tsv`, and
  `combined_protein.tsv` contain cross-sample reports with aligned
  sample spectral-count and intensity columns. `combined_protein.tsv` also
  reports sequence coverage from the union of observed peptides. Aerith
  performs this directly from its in-memory scored PSMs, without external
  workspaces or `.meta` intermediates.
- `combined_protein_with_PSM.tsv` is also written natively by Aerith. For each
  sample it contains aligned comma-separated PSM IDs, peptide sequences, MS1
  and MS2 isotopic abundances, and chromatographically quantified SIP
  intensities. Values at the same comma-separated position belong to the same
  accepted PSM.
  Shared peptides and all of their aligned PSM values are included only in
  the selected razor protein row; alternative compatible proteins are not
  given duplicate peptide or intensity values.
  Accepted match-between-runs features are appended to the acceptor sample
  with `<PeptideSequence>_MBR` as the PSM ID, donor MS1/MS2 isotopic
  abundances, and the transferred intensity quantified in the acceptor run.
  It omits `<sample>_ProteinAbundance` and
  `<sample>_log10_precursorIntensities`. When SIP negative controls are
  requested, Aerith also writes
  `combined_protein_with_SIP_filtered_PSM.tsv`; there is no Python report
  implementation or fallback.
- FASTA search stages pass the original target FASTA and generated decoy FASTA
  to Aerith separately; the workflow does not create `targetDecoy.faa`. The
  fast-SIP spectra-search stage gets decoy protein evidence directly from its
  decoy PIN files and uses only the target FASTA for report annotations, so it
  does not generate `spectra_search/decoy.faa`. Standalone spectra-search mode
  runs Aerith filtering only and produces the per-sample
  `*_filtered_psms.tsv` files without protein assembly outputs.
- Before reranking, Aerith removes decoy PSMs whose stripped peptide sequence
  is also present in any target PSM from the input samples and reports the
  removed collision count in the workflow log.

Aerith writes its timing and result summary to the workflow's existing
`sipros_workflow.log`; it does not create a second `aerith.log` in the output
directory. Every stage reports wall time, CPU time, and observed speedup. The
results section reports distinct naked peptide sequences, distinct modified
peptide forms, PTM peptide forms, and PTM-bearing PSMs at the configured FDR.
If one sample has no confident target PSMs with which to initialize its
sample-specific SVM, Aerith warns, writes zero accepted PSMs for that sample,
and continues filtering the remaining samples.

Aerith converts `log10_precursorIntensities` back to linear PSM intensity and
uses a top-three peptide-ion rule for total, unique, and razor protein
intensity in its standard protein reports.

### Native IonQuant-style chromatographic quantification

For FASTA workflows, Aerith reads the MS1 scans directly from each Raxport
HDF5 file. Regular search retains its monoisotopic, M+1, and M+2 traces. SIP
search instead calculates the source-aware exact isotope distribution for the
modified peptide at its reported MS2 labeling abundance and traces the most
probable `--quant-top-isotopes` peaks (six by default). It resamples the XICs
to a uniform retention-time grid, applies Savitzky-Golay smoothing for peak
detection, and subtracts boundary background from the unsmoothed interpolated
trace. The most probable theoretical isotope trace defines the reported
feature geometry; automatic intensity sums the background-corrected
apex-1/apex/apex+1 values across the traced envelope. Aerith reports the
feature apex, apex scan,
start/end retention times, FWHM, traced scans, and intensity in `psm.tsv`,
`ion.tsv`, and `combined_ion.tsv`. This runs in memory between Aerith filtering
and native protein assembly; it does not create a temporary quantification
table.

Aerith also performs IonQuant-style match between runs. It ranks up to ten
donor runs using shared-ion retention-time and intensity Spearman
correlations, aligns each donor locally with the median retention-time shift
and median absolute deviation, traces target and shifted-decoy isotope
envelopes, scores intensity/isotope-fit/mass-error/retention-time features, and
controls transferred ions by posterior ion FDR. Transferred values are marked
`MBR` in the sample and combined ion, modified-peptide, and peptide tables;
missing values are marked `unmatched`.

The LFQ defaults are 10 ppm MS1 tolerance, a 0.4 minute retention-time window,
at least two isotopes, at least three MS1 scans, and cross-run intensity
normalization. They can be adjusted with
`--quant-mz-ppm`, `--quant-rt-window`, `--quant-min-isotopes`, and
`--quant-min-scans`. `--quant-intensity-mode` accepts IonQuant-compatible
`0` (background-corrected apex), `1` (background-corrected area), or `2`
(IonQuant automatic conventional-LC apex selection, the default). Use
`--no-quant-normalization` to retain raw feature intensities in combined
reports. MBR is controlled with `--mbr-rt-window`, `--mbr-top-runs`,
`--mbr-min-correlation`, and `--mbr-ion-fdr`, or disabled with `--no-mbr`.

`combined_modified_peptide.tsv`, `combined_peptide.tsv`, and
`combined_protein.tsv` contain sample-specific `MaxLFQ Intensity` columns.
Aerith computes median pairwise log ion ratios, solves each connected sample
graph, and anchors its absolute scale to the contributing ion intensities.

The native C++ workflow in `wf33` launches one cross-sample Aerith process for
filtering, chromatographic
quantification, normalization, MaxLFQ, and protein assembly. The optional SIP
negative-control pass selects only primary target PSMs that already passed the
normal target/decoy filter, relabels the designated control samples, and
rescores them inside that same Aerith process using the in-memory dataset.
The secondary SVM reuses each PSM's primary `delta_RT_loess` value and does
not rerun RT prediction or calibration. It does not perform a second
peptide-collision filter, write
or reread `SIP.pin`, launch another Aerith process, or use a Python filtering
filtering. The Aerith report includes the total and per-stage native
negative-control timings, the configured SIP label threshold, counts before
and below that threshold, and the three-fold secondary SVM feature weights.
There is no separate `quant.py` or IonQuant command in the workflow.

### Native predicted-spectrum entropy in Aerith

Aerith can predict DIA-NN fragment spectra directly through LibTorch and add
MSBooster-compatible `unweighted_spectral_entropy` to SVM rescoring. The
renamed DIA-NN 2.6.1 fragmentation checkpoint is
`bin/diann-2.6.1-fragmentation.pt`; the exact DIA-NN token dictionary is
compiled into Aerith.
The implementation keeps DIA-NN's predicted top fragments, applies each
Raxport scan's MS2 m/z window, matches experimental peaks at 20 ppm by default,
and computes the entropy similarity when at least two predicted ions have
experimental matches.
For SIP FASTA PSMs, Aerith expands every predicted b/y product ion through the
same source-aware theoretical isotope model used by Sipros, shifts each peak
according to fragment charge, and weights it by its exact isotope probability
before entropy matching. The selected isotope probabilities are renormalized
per product ion so their expanded intensities sum to the original DIA-NN
product-ion intensity. Aerith retains the five most probable isotopes of every
product ion by default, with no isotope-probability cutoff and no
spectrum-wide 20-peak limit. Change the per-ion count with
`--product-top-isotopes`. All retained DIA-NN product ions are expanded.

Build the release binaries with:

```bash
./make.sh build
```

For `build` and `package` only, the script creates or reuses the
`sipros5-release` micromamba environment. It pins the glibc 2.17 sysroot,
non-MPI HDF5 2.x, PyTorch 2.12.1 CPU/MKL, and Dear ImGui 1.92.9 builds.
Sipros, Aerith, and `siproswf`
dynamically link their external Conda libraries. `package` bundles the complete
non-glibc runtime closure, including HDF5, OpenMP, compiler runtimes, Torch,
MKL, ImGui, GLFW, and OpenGL support libraries, under `bin/lib`.
The resulting binaries require host glibc 2.17 or newer. Neither command
configures, builds, or packages `siprosMPI`.

The optional Conda GPU build remains available from the `sipros5` environment:

```bash
./make.sh buildConda
```

The Aerith Conda build requires PyTorch/LibTorch 2.12.1 or newer. For the GPU build, install the matching CUDA 12.9
PyTorch and CUDA development packages shown at the top of `make.sh`.

`make.sh` obtains `Torch_DIR` from the environment's Python package and checks
that Aerith dynamically resolves `libtorch` and `libc10`; LibTorch is not
statically embedded in Aerith.

Aerith automatically uses CUDA for DIA-NN spectrum and RT prediction when a
CUDA-enabled PyTorch installation and an accessible GPU are available. If CUDA
is unavailable or initialization/inference fails, Aerith retries prediction on
the CPU. The FASTA workflow does not override `CUDA_VISIBLE_DEVICES`,
so it honors GPU visibility assigned by the user, container, or scheduler. The
selected device or fallback reason is written to the workflow log. The timing
table reports the exact DIA-NN model-inference wall and CPU times with `(GPU)`
or `(CPU)` in the row label, in addition to the complete feature-stage time.

Pass one Raxport HDF5 file for every input sample, in the same order as the PIN
pairs:

```bash
build/conda/bin/aerith \
  --target-pin sample_target.pin \
  --decoy-pin sample_decoy.pin \
  --spectra sample.h5 \
  --output-prefix results/sample
```

Use `--spectrum-model` to override the checkpoint file and `--fragment-ppm` to
change the matching tolerance. A build configured with
`-DAERITH_ENABLE_TORCH=OFF` remains usable for filtering without `--spectra`.

The complete FASTA workflow supplies each sample's Raxport HDF5 file
to Aerith automatically, so no extra flag is needed:

```bash
siproswf -i raw -f Ecoli.fasta -t 40 -o regular_output
```

Regular and SIP FASTA searches both receive predicted spectrum and RT
features; SIP applies its theoretical product-ion isotope envelopes.
Search-spectra keeps its existing feature set. The DIA-NN license notice
available for the checkpoint is reproduced in `LICENSE`; review the
version-specific terms before redistributing the checkpoint or a package
containing it.

When this feature is enabled, Aerith's timing report includes `Predict spectra
and compute entropy`. This measures the complete feature-generation stage,
including unique peptide-charge preparation, Torch inference, HDF5 loading,
fragment matching, and entropy calculation.

### Native DIA-NN delta-RT in Aerith

Aerith also loads the renamed DIA-NN 2.6.1 `models/rt.d0.pt` checkpoint from
`bin/diann-2.6.1-retention-time.pt`. It predicts each unique modified peptide
sequence once (RT is charge-independent in this graph), converts the raw model
output to DIA-NN iRT, and fits a separate monotonic robust-LOESS calibration
for each input sample. It computes three MSBooster-compatible values:

- `delta_RT_loess`: absolute predicted-iRT error after mapping observed RT to
  iRT. This is the only DIA-NN RT value used for SVM training and shown in the
  SVM feature-weight report.
- `delta_RT_loess_real`: absolute error in the sample's RT units. This is an
  output-only diagnostic written immediately after `retentiontime`.
- `pred_RT_real_units`: predicted iRT mapped back to the sample's RT units.
  This is also output-only and follows `delta_RT_loess_real`.

The calibration uses up to 5,000 unique precursors ranked only by `WDPscores`:
the best 4,000 globally plus up to 20 additional nonduplicate precursors from
each of 50 equal-width observed-RT bins. Sparse bins are backfilled with the
best remaining WDP candidates. Target/decoy labels, `log10_evalue`, and
`hyperscore` are not consulted.
Both Sipros symbol modifications and bracketed numeric mass notation are
accepted. Use `--rt-model` to override the checkpoint. The timing report's
`Predict RT and compute delta-RT` row includes peptide deduplication, Torch
inference, bandwidth selection, LOESS calibration, inverse mapping, and all
three value calculations. When the DIA-NN RT feature is generated,
Aerith skips its legacy nested chemical RT model and does not add the redundant
`sqrtAbsDeltaRT` feature. Use `--no-predicted-rt` to disable DIA-NN RT and use
the legacy Aerith RT model instead.

## Sipros5 Release Setup Guide

### 1. Create Conda Environment

```bash
conda create -n sipros5 lxml pandas seqkit python=3.12 -c bioconda -c conda-forge
conda activate sipros5
```

### 2. Download Sipros5 Release

```bash
wget https://github.com/xyz1396/sipros5/releases/download/6.0.0/sipros_linux_6.0.0.zip
unzip sipros_linux_6.0.0.zip
chmod +x sipros/bin/* sipros/wf33/extractPro.sh
```

### 3. Example Commands

#### Regular Search

```bash
sipros/bin/siproswf --regular-fasta-search -i raw/Pan_062822_X1iso5.raw -f Ecoli.fasta -o regular_output
```

#### Extract protein sequences identified in Regular search

```bash
sipros/wf33/extractPro.sh Ecoli.fasta regular_output/combined_protein.tsv db.faa
```

#### Label Search

```bash
sipros/bin/siproswf --sip-fasta-search -i raw -f db.faa -e C13 -o sip_output
```

#### Label Search with negative control using unlabeled sample

```bash
sipros/bin/siproswf --sip-fasta-search -i raw -f db.faa -e C13 --negative_control Pan_062822_X1iso5 -o sip2_output
```

### 6. Citations

1. Xiong, Y., Mueller, R. S., Feng, S., Guo, X., & Pan, C. (2024). [Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins](https://doi.org/10.1186/s40168-024-01866-1). *Microbiome*, 12, 148.
2. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). [CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing](https://doi.org/10.1093/bioadv/vbae024). *Bioinformatics Advances*, 4(1), vbae024.
3. Guo, X., Li, Z., Yao, Q., Mueller, R. S., Eng, J. K., Tabb, D. L., Hervey IV, W. J., & Pan, C. (2018). [Sipros Ensemble improves database searching and filtering for complex metaproteomics](https://doi.org/10.1093/bioinformatics/btx601). *Bioinformatics*, 34(5), 795-802.
4. Wang, Y., Ahn, T.-H., Li, Z., & Pan, C. (2013). [Sipros/ProRata: a versatile informatics system for quantitative community proteomics](https://doi.org/10.1093/bioinformatics/btt329). *Bioinformatics*, 29(16), 2064-2065.
5. Pan, C., Kora, G., McDonald, W. H., Tabb, D. L., VerBerkmoes, N. C., Hurst, G. B., Pelletier, D. A., Samatova, N. F., & Hettich, R. L. (2006). [ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation](https://doi.org/10.1021/ac060654b). *Analytical Chemistry*, 78(20), 7121-7131.
