
```bash
conda create -n sklearn scikit-learn pandas matplotlib python=3.10
conda activate sklearn
# update py2 to py3
2to3 --output-dir=script3 -W -n script
```