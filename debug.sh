#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

# ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct1/ft/Pan_062822_X1iso5.FT2'

ft='/scratch/yixiong/benchmark/astral/ftPct1/X3_ID110156_01_OA10034_10302_120823.1.FT2'

echo SiprosArgs="-f $ft -fasta fasta/Decoy.fasta -c SiprosEnsembleConfig.cfg -o regular" >.vscode/.env

echo SiprosArgs="-w mzml -fasta fasta/Decoy.fasta -c SiprosEnsembleConfig.cfg -o regular" >.vscode/.env