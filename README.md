# Sipros5 Setup Guide

## 1. Create Conda Environment

```bash
conda create -n sipros5 lxml pandas python=3.12
conda activate sipros5
```

## 2. Download Sipros5 Release

```bash
wget https://github.com/xyz1396/sipros5/releases/download/5.0/siprosRelease.zip
mkdir sipros
unzip siprosRelease.zip -d sipros
chmod +x sipros/tools/* sipros/siproswf
```

## 3. Download Raw Files

```bash
mkdir raw
# Download raw file with 1% 13C
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD041414/Pan_062822_X1iso5.raw -P raw
# Download raw file with 50% 13C  
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD041414/Pan_052322_X13.raw -P raw
```

## 4. Download E. coli Protein FASTA Sequence

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz
gunzip UP000000625_83333.fasta.gz -c > Ecoli.fasta
```

## 5. Example Commands

### Regular Search

```bash
python sipros/script33/main.py -i raw/Pan_062822_X1iso5.raw -f Ecoli.fasta -t 40 -o regular_output
```

### Label Search

```bash
python sipros/script33/main.py -i raw -f Ecoli.fasta -e C13 -t 40 -o sip_output
```

### Label Search with negative control using unlabeled sample

```bash
python sipros/script33/main.py -i raw -f Ecoli.fasta -e C13 -t 40 --negative_control Pan_062822_X1iso5 -o sip2_output
```

## 6. Citation

1. Xiong, Y., Mueller, R.S., Feng, S., Guo, X. and Pan, C., 2024. Proteomic stable isotope probing with an upgraded Sipros algorithm for improved identification and quantification of isotopically labeled proteins. *Microbiome*, 12.
2. Li, J., Xiong, Y., Feng, S., Pan, C., & Guo, X. (2024). CloudProteoAnalyzer: scalable processing of big data from proteomics using cloud computing. *Bioinformatics Advances*, vbae024.
3. Guo, X., Li, Z., Yao, Q., Mueller, R.S., Eng, J.K., Tabb, D.L., Hervey IV, W.J. and Pan, C., 2018. Sipros ensemble improves database searching and filtering for complex metaproteomics. *Bioinformatics*, 34(5), pp.795-802.
4. Wang, Y., Ahn, T.H., Li, Z. and Pan, C., 2013. Sipros/ProRata: a versatile informatics system for quantitative community proteomics. *Bioinformatics*, 29(16), pp.2064-2065.
5. Pan, C., Kora, G., McDonald, W.H., Tabb, D.L., VerBerkmoes, N.C., Hurst, G.B., Pelletier, D.A., Samatova, N.F. and Hettich, R.L., 2006. ProRata: a quantitative proteomics program for accurate protein abundance ratio estimation with confidence interval evaluation. *Analytical Chemistry*, 78(20), pp.7121-