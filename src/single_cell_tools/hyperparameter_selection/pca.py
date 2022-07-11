import scipy as sp
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD

def pca_explained_variance(adata,threshold=0.9):
    """
    Function that chooses the number of PCA components based on the explained variance of the first PCAs.

    **Argumets:**\n

     - **adata**: Annotated Data object with the pca already computed and added the metainformation to .uns.\n
     - **threshold=0.9**: Threshold specifying the amount of explained variance of the retained PCs.\n

    **Returns:**\n

    Adds to adata.obs["pca"] information of the hyperparameter selection method "hps_method" 
    and the number of retained PCs under "n_pca_relevant".\n
    """
    
    aux = np.where(adata.uns["pca"]["variance_ratio"].cumsum() >= threshold)[0]
    if len(aux) == 0:
        adata.uns["pca"]["n_pca_relevant"] = len(adata.uns["pca"]["variance_ratio"].cumsum())
    else:
        adata.uns["pca"]["n_pca_relevant"] = aux[0] - 1

    adata.uns["pca"]["hps_method"] = {"name":"explained_variance"}
    
    return

def pca_permutation_test(adata,p_value=0.05,n_perms=10,key_obsm="X",use_HVGs=False,key_HVGs="highly_variable",strict=True):
    """
    Function that chooses the number of PCA components based on a permutation test.

    **Argumets:**\n

     - **adata**: Annotated Data object with the pca already computed and added the metainformation to .uns.\n
     - **p_value=0.05**: P value retaining the number of PCs under the null model.\n
     - **n_perms=100**: Number of times that the randomization is performed to compute the p-values.\n
     - **key_obsm="X"**: Matrix to be used for computing the permutation. If "X", uses adata.X.\n
     - **use_HVGs=False**: Whether to remove HVGs for computing the randomized PCs. It should be consistent with the computation of the PCs without randomization.\n
     - **key_HVGs="highly_variable"**: Key in .var that contains the information of the Highly varying genes.\n
     - **strict=True**: If True, makes the p-value comparing the most informative random PC. If False, compares each PC component with the random component in the same position. 

    **Returns:**\n

    Adds to adata.obs["pca"] information of the hyperparameter selection method "hps_method",
    the number of retained PCs under "n_pca_relevant", the variance randomized mean over the permutations "variance_randomized_mean"
    and the p values per component under "permutation_p_value".\n
    """
    
    if key_obsm == "X":
        X = adata.X
    else:
        X = adata.obsm["X"]
        
    if use_HVGs:
        X = X[:,adata.var[key_HVGs]]
    X = X.copy()
        
    if not "pca" in adata.uns.keys():
        ValueError("A pca algorithm has to be runned before.")
        
    length = adata.uns["pca"]["variance"].shape[0]
    if sp.sparse.issparse(X):
        model = TruncatedSVD(n_components=length)
    else:
        model = PCA(n_components=length)
        
    shuffle_matrix = np.repeat(np.array(range(X.shape[0])).reshape(-1,1),X.shape[1],axis=1)
    variance = np.zeros([n_perms,length])
    for perm in range(n_perms):
        shuffle_matrix = np.apply_along_axis(np.random.permutation,0,shuffle_matrix)
        X_shuffled = X[shuffle_matrix,range(X.shape[1])]
        
        model.fit(X_shuffled)
        
        variance[perm,:] = model.explained_variance_
        
    adata.uns["pca"]["variance_randomized_mean"] = np.mean(variance, axis=0)
    if strict:
        var = np.repeat(np.max(variance[:,:adata.uns["pca"]["variance"].shape[0]],axis=1).reshape(-1,1),length,axis=1)
        adata.uns["pca"]["permutation_p_value"] = np.mean((var > np.repeat(adata.uns["pca"]["variance"].reshape(1,-1),n_perms,axis=0)), axis=0)        
    else:
        adata.uns["pca"]["permutation_p_value"] = np.mean((variance[:,:adata.uns["pca"]["variance"].shape[0]] > adata.uns["pca"]["variance"].reshape(1,-1)), axis=0)
    pos = np.where(adata.uns["pca"]["permutation_p_value"]>p_value)[0]
    if len(pos) != 0:
        adata.uns["pca"]["n_pca_relevant"] = pos[0]
    else:
        print("Not found a cut point in the first ",length," PCs.")
        adata.uns["pca"]["n_pca_relevant"] = length
        
    adata.uns["pca"]["hps_method"] = {"name":"permutation_test","n_perm":n_perms,"p_value":p_value}
    
    return