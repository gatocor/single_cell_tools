import pandas as pd
import scipy as sp
import numpy as np

def filter_cells(adata,*args):
    """
    Function to filter cells by obs or var features.
    
    **Arguments:**
     - **adata**: Annotated Dataset to filter.
     - **args**: Arguments seting the filters. Can be defined in three ways:
        - Filtering obs categorical features: a tuple of the form (obs,[retained_category1,retaiden_category2...])
        - Filtering obs numerical features: a tuple of the form (obs,min,max)
        - Filtering var features: a tuple of the form (var_key,gene,min,max)
        
    **Returns:**
     - Array of length of cells with True in the retained features.
     - DataFrame with statistics of removed objects.
    """
    
    data = pd.DataFrame(columns=["Object","#below","%below","#above","%above","Total_removed","%Total_removed"])
    N = adata.obs.shape[0]
    
    retained = np.array([True for i in range(N)])
    for i in range(len(args)):
        filt = args[i]
        if len(filt) == 2:
            retained_categories = [i in filt[1] for i in adata.obs[filt[0]]]
            
            retained *= retained_categories
            
            data.loc[i,:] = ["obs."+filt[0],
                         0,0,
                         0,0,
                         np.sum(np.invert(retained_categories)),np.sum(np.invert(retained_categories))/N]
        elif len(filt) == 3:
            upper = adata.obs[filt[0]] >= filt[1]
            lower = adata.obs[filt[0]] <= filt[2]
            
            retained *= upper
            retained *= lower
            
            data.loc[i,:] = ["obs."+filt[0],
                         np.sum(np.invert(upper)),np.sum(np.invert(upper))/N,
                         np.sum(np.invert(lower)),np.sum(np.invert(lower))/N,
                         np.sum(np.invert(lower*upper)),np.sum(np.invert(lower*upper))/N]
        elif len(filt) == 4:
            gen = np.where(adata.var[filt[0]] == filt[1])[0][0]
            
            if sp.sparse.issparse(adata.X):
                upper = np.array((adata.X[:,gen].todense() >= filt[2]))[:,0]
                lower = np.array((adata.X[:,gen].todense() <= filt[3]))[:,0]
            else:
                upper = np.array((adata.X[:,gen] >= filt[2]))[:,0]
                lower = np.array((adata.X[:,gen] <= filt[3]))[:,0]
            
            retained *= upper
            retained *= lower
            
            data.loc[i,:] = ["var."+filt[0]+".X."+filt[0],
                         np.sum(np.invert(upper)),np.sum(np.invert(upper))/N,
                         np.sum(np.invert(lower)),np.sum(np.invert(lower))/N,
                         np.sum(np.invert(lower*upper)),np.sum(np.invert(lower*upper))/N]

    data.loc[len(args),:] = ["TOTAL_",0,0,0,0,np.sum(np.invert(retained)),np.sum(np.invert(retained))/N]
        
    return retained.values, data