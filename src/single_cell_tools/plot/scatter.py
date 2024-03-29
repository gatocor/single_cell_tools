import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches

def plot_scatter_labels(ax,adata,obsm,obs,**kwargs):

    X = adata.obsm[obsm]
    for label in np.unique(adata.obs[obs].values):

        subset = adata.obs[obs].values == label
        pos = np.mean(X[subset],axis=0)
        ax.text(pos[0],pos[1],label,**kwargs)

    return

def plot_scatter_lines(ax,adata,obsm,obs,matrix,width=10):

    categories = adata.obs[obs]
    for i,origin in enumerate(categories.cat.categories.values): 
        for j,target in enumerate(categories.cat.categories.values):
            x1,y1=adata[categories==origin].obsm[obsm].mean(axis=0)
            x2,y2=adata[categories==target].obsm[obsm].mean(axis=0)
            
            ax.plot([x1,x2],[y1,y2],linewidth=width*matrix[i,j],color="black")

    return

def plot_scatter_pies(ax,adata,obsm,obs,obs2,width=10,cmap="rocket_r"):

    categories = adata.obs[obs]
    categories2 = np.sort(adata.obs[obs2].cat.categories.values)
    color = {i:j for i,j in zip(np.sort(categories2),sns.color_palette(cmap,len(categories2)))}
    for i,origin in enumerate(categories.cat.categories.values):
            xmin,ymin = adata.obsm[obsm].min(axis=0)
            xmax,ymax = adata.obsm[obsm].max(axis=0)
            x1,y1=adata[categories==origin].obsm[obsm].mean(axis=0)
            subax = ax.inset_axes([(x1-xmin)/(xmax-xmin)/1.1,(y1-ymin)/(ymax-ymin)/1.1,.1,.1])
            m = adata[categories==origin].obs.groupby(obs2).count().iloc[:,0]
            cat = m.index
            c = [color[i] for i in m.index]
            subax.pie(m.values,colors=c)

    patches = [mpatches.Patch(color=color[k], label=k) for k in color.keys()]
    ax.legend(title=obs2,fontsize=20,title_fontsize=20,handles=patches)

    return