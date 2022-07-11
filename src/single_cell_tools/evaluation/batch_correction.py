from sklearn.neighbors import NearestNeighbors
import numpy as np
import scipy as sp

def kBET(rep, labels, pReject = 0.01, nTopVectors=None, knns = [10,20,30], *args, **kwargs):
    """
    Function to compute the kBET metric for batch correction mixing.
    
    **Argumnets:** 
     - **rep**: Representation to use for computing the KNN graph.
     - **labels**: Labels specifying the batch of each cell.
     - **pReject = 0.01**: Rejection rate for a cell to have a mixed neighborhood.
     - **nTopVectors=None**: Number of first columns to use from the rep to compute the KNN. if `None`, use all of them.
     - **knns=[10,20,30]**: List of knnn to be used.
     - **args**: Args to be send to sklearn.NearestNeighbors function.
     - **kwargs**: Keyword Args to be send to sklearn.NearestNeighbors function.
        
    **Return:**
     Median over the mean rejection rate of each knns.
     P value for each samples, it helps to find potential clusters that did not mixed well.
    """
    
    if nTopVectors != None:
        repEffective = rep[:,:nTopVectors]
    else:
        repEffective = rep
    
    size = rep.shape[0]
    mpj = np.zeros(len(knns))
    for j,knn in enumerate(knns):
        #Fit the model
        model = NearestNeighbors(knn,*args,**kwargs).fit(rep)
        #Compute neighbors
        m = model.kneighbors()[1]
        #Compùte the nij and fi x k of the paper to compùte the X^2 distribution
        l = np.unique(labels)
        nij = np.zeros([len(labels),len(l)])
        n = np.zeros(len(l))
        for i,sample in enumerate(l):
            nij[:,i] = np.sum(labels[m]==sample,axis=1)
            n[i] = np.sum(labels == sample)/len(labels)*knn

        #Compute the kij
        kj = np.sum((nij-n)**2/n,axis=1)
        #Compute the cumulative
        pj = 1-sp.stats.chi2(len(l)).cdf(kj)
        #Make the mean between all subsamples
        mpj[j] = np.mean(pj > pReject)
    
    return np.median(mpj), pj