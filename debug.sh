#!/bin/bash

ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct1/ft/Pan_062822_X1iso5.FT2'

# ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/sipros5/test/mzml/Pan_062822_X1iso5.mzML'

# ft='/scratch/yixiong/benchmark/astral/ftPct1/X3_ID110156_01_OA10034_10302_120823.1.FT2'

# ft='/scratch/yixiong/Astral/K14/K14_ID116311_01_SPD30_OA10034_10537_071524.100208.FT2'

# echo SiprosArgs="-f $ft -fasta fasta/Decoy.fasta -c SiprosEnsembleConfig.cfg -o tmp" >.vscode/.env
# echo SiprosArgs="-f $ft -fasta fasta/Decoy.fasta -c wf_output/configs/SIP_C13_02.000Pct.cfg -o tmp" >.vscode/.env

# echo SiprosArgs="-w mzml -fasta fasta/Decoy.fasta -c SiprosEnsembleConfig.cfg -o regular" >.vscode/.env

cat <<EOL >.vscode/.env
SiprosArgs="-f $ft -fasta fasta/Decoy.fasta -c wf_output/configs/SIP_C13_02.000Pct.cfg -o tmp"
AerithArgs="-t wf_outputAstral50/K15_ID116312_01_SPD30_OA10034_10537_071524/target \
-d wf_outputAstral50/K15_ID116312_01_SPD30_OA10034_10537_071524/decoy \
-n 8 -f wf_outputAstral50/K15_ID116312_01_SPD30_OA10034_10537_071524/ft \
-j 10 -c /nullspace/sipros5/configTemplates/SIP.cfg \
-p 0 -o wf_outputAstral50/K15_ID116312_01_SPD30_OA10034_10537_071524/K15_ID116312_01_SPD30_OA10034_10537_071524.pin"
OMP_STACKSIZE=16M
EOL