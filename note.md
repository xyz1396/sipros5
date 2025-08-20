
# update py2 to py3

```bash
conda create -n sklearn scikit-learn pandas matplotlib python=3.10
conda activate sklearn
# update py2 to py3
2to3 --output-dir=script3 -W -n script
```

# run test

```bash
setsid ./test.sh runSIP > test/test.log 2>&1 &
setsid ./test.sh runSIPone > test/test.log 2>&1 &
setsid ./test.sh run > test/test.log 2>&1 &
```

# convert mzml

```bash
mono test/tools/ThermoRawFileParser/ThermoRawFileParser.exe -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/raw/pct1/Pan_062822_X2iso5.raw -o test/mzml
```
# convert ft

```bash
test/tools/Raxport -f /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/raw/pct1/Pan_062822_X2iso5.raw -o test/ft
```

# make full decoy

```bash
python script3/makeReverseDecoyDB.py test/fasta/EcoliWithCrapNodup.fasta test/fasta/FullDecoy.fasta
```

# test AerithFeatureExtractor 
```bash
setsid python script33/debugProxy.py > test/test.log 2>&1 &

bin/AerithFeatureExtractor \
-t /nullspace/sipros5/test/wf_output/Pan_052322_X13/target \
-d /nullspace/sipros5/test/wf_output/Pan_052322_X13/decoy \
-n 8 \
-f /nullspace/sipros5/test/wf_output/Pan_052322_X13/ft \
-j 8 \
-c /nullspace/sipros5/configTemplates/SIP.cfg \
-p 0 \
-o /nullspace/sipros5/test/wf_output/Pan_052322_X13/X13.pin
```

# test sipros_wf

```bash
cd /scratch/yixiong/recovered/benchmark/cecum_13C_glucose
/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/sipros5/siproswf \
        -i raw \
        -f db.faa \
        -e C13 \
        -t 40 \
        -n 6 \
        -o cecum_wf_SIP \
        --dryrun \
        --negative_control "M1_PanC_20250327_01_DDA_01,M2_PanC_20250327_01_DDA_02,M3_PanC_20250327_01_DDA_03" \
        > cecum.log 2>&1 

rm -r test/wf_output
mkdir test/wf_output

setsid ./siproswf \
        -i /prebiotics/ubuntuShare/EcoliSIP/goodResults/pct1/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -t 160 \
        -o test/wf_outputR \
> testR.log 2>&1 &


cd /ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/sipros5

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct1/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        --dryrun \
        -o test/wf_output1 \
> test/test1.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct2/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        --dryrun \
        -o test/wf_output2 \
> test/test2.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct5/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        --dryrun \
        -o test/wf_output5 \
> test/test5.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct25/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        -o test/wf_output25 \
> test/test25.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct50/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        --dryrun \
        -o test/wf_output50 \
> test/test50.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct99/raw \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        --dryrun \
        -o test/wf_output99 \
> test/test99.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/AstralC13/pct1 \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -t 40 \
        -n 15 \
        -s \
        -o test/wf_AstralR1 \
> test/testAstralR1.log 2>&1 &

setsid ./siproswf \
        -i /ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/AstralC13/pct50 \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 160 \
        -s \
        -n 15 \
        -o test/wf_outputAstral50 \
> test/testAstral50.log 2>&1 &

./siproswf \
        -i /ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/AstralC13/pct99 \
        -f test/fasta/EcoliWithCrapNodup.fasta \
        -e C13 \
        -t 40 \
        -s \
        -n 15 \
        -o test/wf_outputAstral99 \
        --dryrun
```
