import numpy as np

def downsample(adata,p,batch=None,replace=False,seed=None):
    """
    Function to downsample an Annotated data. If a categorical batch label is givven the downsample is performed over each label independently.
    This can be useful when for some reason you want to keep small clusters present in the data.

    **Arguments:**
     - **adata**: Annotated Data to downsample.
     - **p**: Proportion of kept cells.
     - **batch=None**: Categorical label in .obs to use for downsampling by category.
     - **replace=False**: If to downsample with replacement or not.
     - **seed=None**: Random seed to start the downsampling.

    **Returns:**
     Downsampled copy of the Annotated Data
    """
    
    if seed != None:
        np.random.seed(seed)
    
    if batch == None:
        N = int(np.ceil(adata.shape[0]*p))
        l = np.random.choice(adata.obs.index,N,replace=replace)
    else:
        l = []
        for category in adata.obs[batch].cat.categories.values:
            N = int(np.ceil(adata.obs[adata.obs[batch]==category].shape[0]*p))
            l = np.append(l,np.random.choice(adata.obs[adata.obs[batch]==category].index,N,replace=False))
        
        l = adata.obs.index.isin(l)
        
    return adata[l,:]