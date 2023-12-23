#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct1/ft/Pan_062822_X1iso5.FT2'

echo SiprosArgs="-f $ft -fasta fasta/FullDecoy.fasta \
    -c SiprosEnsembleConfig.cfg -o regular" >.vscode/.env
