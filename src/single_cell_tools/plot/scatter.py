import numpy as np
import seaborn as sns

def plot_scatter(adata,ax,obsm,hue=None,reverse=False,**kwargs):

    if hue != None:
        color = adata.obs[hue]
        X = adata.obsm[obsm]

        order = np.argsort(color)

        sns.scatterplot(x=X[order,0],y=X[order,1],hue=color,**kwargs)
    
    else:
    
        X = adata.obsm[obsm]

        sns.scatterplot(x=X[:,0],y=X[:,1],**kwargs)

    return

def plot_scatter_labels(adata,ax,obsm,obs,**kwargs):

    X = adata.obsm[obsm]
    for label in adata.obs[obs].unique().values:

        subset = adata.obs[obs] = label
        pos = np.mean(X[subset])
        ax.text(pos[0],pos[1],label,**kwargs)

    return