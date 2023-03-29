#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd
import scanpy as sc
import sklearn as sk
from anndata import AnnData
from numbers import Number
import warnings
from typing import Union, Optional, Tuple, Collection, Sequence, Iterable
import scipy
from scipy.spatial import distance
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from uclab import uclab, uclab_split
import pyarrow as pa
from pyarrow import ChunkedArray
from sk.decomposition import TruncatedSVD

def scsampler(
    data: Union[AnnData, np.ndarray, spmatrix,ChunkedArray],
    fraction: Optional[float] = None,
    n_obs: Optional[int] = None,
    random_state: int = 0,
    copy: bool = False,
    obsm: Optional[str] = 'X_pca',
    dr_num: Optional[int]=None,
    obs_index: int = 0,
    var_index: int = 0,
    random_split: Optional[int] = None,
) -> Optional[AnnData]:
    """    Subsample to a fraction of the number of observations. This function refers to the subsample function in scanpy. 
    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    fraction
        Subsample to this `fraction` of the number of observations.
    n_obs
        Subsample to this number of observations.
    random_state
        Random seed to change subsampling.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.
    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~anndata.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """

    """    Subsample to a fraction of the number of observations.  
    Parameters
    ----------
    data
        pyarrow.ChunkedArray
    fraction
        Subsample to this `fraction` of the number of observations.
    n_obs
        Subsample to this number of observations.
    random_state
        Random seed to change subsampling.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.
    dr_num
       If data format is pyarrow.ChunkedArray,
       The number of dimensions.
    obs_index
       Reset the cell indeces to go from 0 to N (number of cells in the query) 
    var_index
       Reset the gene indeces to go from 0 to N (number of genes in the query)
    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~anndata.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """              
    ## Sample space
    if isinstance(data, AnnData):
        X = data.obsm[obsm]
        old_n_obs = data.n_obs
        old_n_vars = X.shape[1]
    elif isinstance(data, ChunkedArray):
        X = data.to_numpy()
        X = scipy.sparse.coo_matrix((X, (obs_index, var_index))).tocsr()
        old_n_obs = X.shape[0]
        old_n_vars = X.shape[1]
    else: 
        data
        old_n_obs = data.shape[0]
        old_n_vars = X.shape[1]
    
    np.random.seed(random_state)
    
    if dr_num is not None:
        svd = TruncatedSVD(dr_num)
        X = svd.fit_transform(X)
    if n_obs is not None:
        new_n_obs = n_obs
    elif fraction is not None:
        if fraction > 1 or fraction < 0:
            raise ValueError(f'`fraction` needs to be within [0, 1], not {fraction}')
        new_n_obs = int(fraction * old_n_obs)
    else:
        raise ValueError('Either pass `n_obs` or `fraction`.')
    
    ## Run uclab
    split = 1 if random_split is None else random_split
    obs_indices = uclab_split(X, new_n_obs, alpha = 4 * old_n_vars, drop_start=1, drop_rate=0, split=split)
    
    if isinstance(data, AnnData):
        if copy:
            return data[obs_indices].copy()
        else:
            data._inplace_subset_obs(obs_indices)
    elif isinstance(data, ChunkedArray):
        if copy:
            return X[obs_indices,:].copy(),obs_indices
        else:
            return X[obs_indices,:],obs_indices
    else:
        X = data
        if copy:
            return X[obs_indices], obs_indices
        else:
            return obs_indices


# In[ ]:




