import scipy as sp
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD

def explained_variance_pca(adata,threshold=0.9,key_obsm="X"):
    
    aux = np.where(adata.uns["pca"]["variance_ratio"].cumsum() >= threshold)[0]
    if len(aux) == 0:
        adata.uns["pca"]["n_pca_relevant"] = len(adata.uns["pca"]["variance_ratio"].cumsum())
    else:
        adata.uns["pca"]["n_pca_relevant"] = aux[0] - 1
        
    adata.uns["pca"]["hps_method"] = {"name":"explained_variance"}
    
    return

def permutation_test_pca(adata,p_value=0.05,n_perms=100,key_obsm="X",use_HVGs=True,key_HVGs="highly_variable"):
    
    if key_obsm == "X":
        X = adata.X
    else:
        X = adata.obsm["X"]
        
    if use_HVGs:
        X = X[:,adata.obs[key_HVGs]]
    X = X.copy()
        
    if not "pca" in adata.uns.keys():
        ValueError("A pca algorithm has to be runned before.")
        
    if sp.sparse.issparse(X):
        model = TruncatedSVD()
    else:
        model = PCA()
        
    shuffle_matrix = np.repeat(np.array(range(X.shape[0])).reshape(-1,1),X.shape[1],axis=1)
    variance = np.zeros([n_perms,X.shape[1]])
    for perm in range(n_perms):
        shuffle_matrix = np.apply_along_axis(np.random.permutation,0,shuffle_matrix)
        X_shuffled = X[shuffle_matrix,range(X.shape[1])]
        
        model.fit(X_shuffled)
        
        variance[perm,:] = model.explained_variance_
        
    adata.uns["pca"]["variance_randomized_mean"] = np.mean(variance, axis=0)
    adata.uns["pca"]["permutation_p_value"] = np.mean((variance[:,:adata.uns["pca"]["variance"].shape[0]] > adata.uns["pca"]["variance"].reshape(1,-1)), axis=0)
    adata.uns["pca"]["n_pca_relevant"] = np.where(adata.uns["pca"]["permutation_p_value"]>p_value)[0][0]
    adata.uns["pca"]["hps_method"] = {"name":"permutation_test","n_perm":n_perms,"p_value":p_value}
    
    return