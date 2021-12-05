import numpy as np
from scipy.spatial import distance
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
    
def uclab(X, n, alpha, drop_start=1, drop_rate=0):
    N = X.shape[0]
    sample_index = np.full(n, -1, dtype=int)
    distances = np.full(N, -1)
    
    # Define function
    #int64 = np.int64
    
    # step 0
    initial_index = np.random.randint(N, size=1)
    sample_index[0] = initial_index
    distances = distance.cdist(np.expand_dims(X[sample_index[0],:], axis=0), X).flatten()
    distances[initial_index] = -1
    d_index = np.argwhere(distances != -1).flatten()
    distances[d_index] = np.power(1/distances[d_index], alpha)
    
    for i in np.arange(1,n):
        #print(i, end=' ')
        # drop large points
        if i == np.int64(drop_start*n) and drop_rate != 0:
            current_number = np.count_nonzero(distances!=-1)
            drop_number = np.int64(current_number*drop_rate)        
            left_number = current_number - drop_number
            more_number = n - i
            # if drop too many, then return current best
            if more_number > left_number:
                drop_number = current_number - more_number
                sort_index = np.argsort(distances) 
                drop_index = sort_index[-drop_number:]
                distances[drop_index] = -1
                keep_index = np.where(distances != -1)[0]
                sample_index[i:] = keep_index
                return sample_index      
            sort_index = np.argsort(distances)            
            drop_index = sort_index[-drop_number:]
            distances[drop_index] = -1
        closest_index = np.where(distances > 0, distances, np.inf).argmin()
        sample_index[i] = closest_index
        distances[closest_index] = -1
        d_index = np.argwhere(distances != -1).flatten()
        d = distance.cdist(np.expand_dims(X[sample_index[i],:], axis=0), X[d_index,:]).flatten()
        d = np.power(1/d, alpha)
        distances[d_index] += d     
    return sample_index


def uclab_split(X, n, alpha, drop_start=1, drop_rate=0, split=4):
    # random shuffle index
    #np.random.seed(seed)
    index = np.arange(X.shape[0])
    np.random.shuffle(index)
    X = X[index,]
    # subsample per split
    X_list = np.array_split(X, split, axis=0)
    n_list = np.array_split(np.arange(n), split, axis=0)
    sample_index_list = [uclab(d, len(n_split), alpha, drop_start, drop_rate) for n_split, d in zip(n_list, X_list)]    
    # adjust split index
    sample_index_list_adjust = []
    adjust = 0
    for i, l in enumerate(sample_index_list):
        if i > 0:
            adjust += len(X_list[i-1])
            l = l + adjust
        sample_index_list_adjust.append(l)
    sample_index = np.concatenate(sample_index_list_adjust)
    # get original index
    sample_index = index[sample_index]
    return sample_index