
# update py2 to py3

```bash
conda create -n sklearn scikit-learn pandas matplotlib python=3.10
conda activate sklearn
# update py2 to py3
2to3 --output-dir=script3 -W -n script
```

# run test

```bash
nohup ./test.sh runSIP > test/test.log 2>&1 &
nohup ./test.sh runSIPone > test/test.log 2>&1 &
nohup ./test.sh run > test/test.log 2>&1 &
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
nohup python script33/debugProxy.py > test/test.log 2>&1 &

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
