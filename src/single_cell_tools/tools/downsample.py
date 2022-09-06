import numpy as np

def downsample(adata,p,batch=None,seed=0):
    
    np.random.seed(seed)
    
    if batch == None:
        N = int(np.ceil(adata.shape[0]*p))
        l = np.random.choice(adata.obs.index,N,replace=False)
    else:
        l = []
        for category in adata.obs[batch].cat.categories.values:
            N = int(np.ceil(adata.obs[adata.obs[batch]==category].shape[0]*p))
            l = np.append(l,np.random.choice(adata.obs[adata.obs[batch]==category].index,N,replace=False))
        
        l = adata.obs.index.isin(l)
        
    return adata[l,:]