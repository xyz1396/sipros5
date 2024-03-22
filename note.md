
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
