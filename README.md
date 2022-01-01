# scSampler

## Overview
`scSampler` is a Python pacakge for fast diversity-preserving subsampling of large-scale single-cell transcriptomic data.

## Installation
Please install it from PyPI:
``` python
pip install scsampler
```

## Quick start
First we load all modules.
```python
import numpy as np
import pandas as pd
import scanpy as sc
from time import time
from scsampler import scsampler
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
```
### Read in data
The example data can be downloaded from https://doi.org/10.5281/zenodo.5811787 in the `anndata` format by `scanpy`. Here we use the ~68'000 PBMC cells. Please modify the path as your own path.
```{python}
adata = sc.read_h5ad('/home/dongyuan/scSampler/data/final_h5ad/pbmc68k.h5ad')
```

### anndata as input
Subsample 10% cells and return a new anndata. The space is top PCs.
```{python}
adata_sub = scsampler(adata, fraction = 0.1, copy = True) 
```
If you want to speed it up, you can use the `random_split`. It will lead to slightly less optimal result, of course.
```{python}
start = time()
adata_sub = scsampler(adata, fraction = 0.1, obsm = 'X_pca', copy = True, random_split = 16)
end = time()
print(end - start)
```
### matrix as input
You can also use the `numpy.ndarray` as the input.
```{python}
mat = adata.obsm['X_pca']
print(type(mat))
res = scsampler(mat, fraction = 0.1, copy = True, random_split = 16)
subsample_index = res[1]
subsample_mat = res[0]
```

## Contact
Any questions or suggestions on `scSampler` are welcomed! If you have any questions, please report it on [issues](https://github.com/SONGDONGYUAN1994/scsampler/issues) or contact Dongyuan (<dongyuansong@ucla.edu>).
