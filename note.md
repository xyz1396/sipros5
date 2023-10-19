
# update py2 to py3

```bash
conda create -n sklearn scikit-learn pandas matplotlib python=3.10
conda activate sklearn
# update py2 to py3
2to3 --output-dir=script3 -W -n script
```

```bash
nohup ./test.sh runSIP > test/test.log 2>&1 &
nohup ./test.sh run > test/test.log 2>&1 &
```
