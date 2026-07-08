#!/usr/bin/env python3
"""Legacy pct benchmark entrypoint.

The FT text benchmark workflow has been removed. Use the universal HDF5
search-spectra workflow through script33/main.py, for example:

  python script33/main.py -i data/pct1/raw/Pan_062822_X1iso5.raw \
      -f data/EcoliWithCrapNodup.fasta -e C13 --range 0-5 --precision 1 \
      --psm-tsv data/pct1/frag/Pan_062822_X1iso5/psm.tsv \
      --unlabeled-input data/pct1/raw/Pan_062822_X1iso5.raw \
      -o data/tmp/raxport_hdf5_workflow_test/search_spectra
"""

from __future__ import annotations

raise SystemExit(__doc__)
