<p align="center">
  <img src="script33/sipros_logo.png" alt="Sipros5 logo" width="240">
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
siproswf -i raw/Pan_062822_X1iso5.raw -f Ecoli.fasta -o regular_output
```

#### Extract protein sequences identified in Regular search

- This step is particularly useful when your protein FASTA is large (for example, several GB in metaproteomics studies). The `regular_output/protein.tsv` file can be replaced with results from other proteomics search engines (e.g., FragPipe, MaxQuant, or Proteome Discoverer) as long as the first column contains the protein identifier.
- If you are working with a small FASTA, you can skip this extraction step and use the original FASTA for the label search.

```bash
extractPro Ecoli.fasta regular_output/protein.tsv db.faa
```

#### Label Search

```bash
siproswf -i raw -f db.faa -e C13 -o sip_output
```

#### Label Search with negative control using unlabeled sample

```bash
siproswf -i raw -f db.faa -e C13 --negative_control Pan_062822_X1iso5 -o sip2_output
```

### 5. Output Files

- `SIP_filtered_psms.tsv`: PSMs from all samples that pass the unlabeled negative-control filter (1% FDR), with SIP element labeling percentages (`MS1IsotopicAbundances`, `MS2IsotopicAbundances`). `isotopicPeakNumbers` is the raw number of extracted MS1 isotope peaks. `MS1IsotopeFitScore` is the theoretical-envelope coverage from `0` to `1`, set to `0` unless at least two compatible peaks are present; scores of at least `0.02` pass the MS1 abundance-fit validity threshold. MS1IsotopicAbundances are more sensitive; MS2IsotopicAbundances are more accurate.
- `protein_with_SIP_filtered_PSM.tsv`: maps unlabeled negative-control filtered PSMs to the proteins identified in each sample.
- For each raw-file subdirectory:
  - `<sample>.h5`: Raxport scan data.
  - `<sample>_target.pin`, `<sample>_decoy.pin`: target and decoy search intermediates.
  - `<sample>.pin`: merged target/decoy PIN consumed by Percolator.
  - `peptide.tsv`, `protein.tsv`: peptides, and proteins passing 1% FDR decoy filtering.
  - `*_filtered_psms.tsv`: PSMs passing 1% FDR decoy filtering with `isotopicPeakNumbers`, `MS1IsotopeFitScore`, `MS1IsotopicAbundances`, and `MS2IsotopicAbundances`.

## Sipros5 Setup Guide (set the python and binary by yourself)

### 1. Create Conda Environment

```bash
conda create -n sipros5 lxml pandas seqkit python=3.12 -c bioconda -c conda-forge
conda activate sipros5
```

### 2. Download Sipros5 Release

```bash
wget https://github.com/xyz1396/sipros5/releases/download/5.0.1/siprosRelease.zip
unzip siprosRelease.zip
chmod +x sipros/tools/* sipros/script33/extractPro.sh
```

### 3. Example Commands

#### Regular Search

```bash
python sipros/script33/main.py -i raw/Pan_062822_X1iso5.raw -f Ecoli.fasta -o regular_output
```

#### Extract protein sequences identified in Regular search

```bash
sipros/script33/extractPro.sh Ecoli.fasta regular_output/protein.tsv db.faa
```

#### Label Search

```bash
python sipros/script33/main.py -i raw -f db.faa -e C13 -o sip_output
```

#### Label Search with negative control using unlabeled sample

```bash
python sipros/script33/main.py -i raw -f db.faa -e C13 --negative_control Pan_062822_X1iso5 -o sip2_output
```

### 6. Citations

1. Xiong, Y., Mueller, R. S., Feng, S., Guo, X., & Pan, C. (2024). [Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins](https://doi.org/10.1186/s40168-024-01866-1). *Microbiome*, 12, 148.
2. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). [CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing](https://doi.org/10.1093/bioadv/vbae024). *Bioinformatics Advances*, 4(1), vbae024.
3. Guo, X., Li, Z., Yao, Q., Mueller, R. S., Eng, J. K., Tabb, D. L., Hervey IV, W. J., & Pan, C. (2018). [Sipros Ensemble improves database searching and filtering for complex metaproteomics](https://doi.org/10.1093/bioinformatics/btx601). *Bioinformatics*, 34(5), 795-802.
4. Wang, Y., Ahn, T.-H., Li, Z., & Pan, C. (2013). [Sipros/ProRata: a versatile informatics system for quantitative community proteomics](https://doi.org/10.1093/bioinformatics/btt329). *Bioinformatics*, 29(16), 2064-2065.
5. Pan, C., Kora, G., McDonald, W. H., Tabb, D. L., VerBerkmoes, N. C., Hurst, G. B., Pelletier, D. A., Samatova, N. F., & Hettich, R. L. (2006). [ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation](https://doi.org/10.1021/ac060654b). *Analytical Chemistry*, 78(20), 7121-7131.
