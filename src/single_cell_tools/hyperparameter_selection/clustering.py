import numpy as np
from cdlib import evaluation, NodeClustering
import igraph
import scanpy as scp

def leiden(adata,
    resolutions=np.arange(.1,1.1,.1),
    evaluation_metrics=[evaluation.surprise,evaluation.modularity_density],
    **kwargs):
    """
    Function that constructs over `scp.tl.leiden <https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.leiden.html#scanpy-tl-leiden>`_
    to automatically select the resolution parameter. The optimization algorithm evaluates several metrics for clustering and chooses the best
    resolution among all the metrics using a Friedman Ranking test. The optimization algorithms uses the implementation from 
    `CDLIB <https://cdlib.readthedocs.io/en/latest/index.html>`_.

    **Arguments:**

     - **adata**: Annotated Data object to cluster. It requires "neighbors" to be already computed and present in .uns.
     - **resolutions=np.linspace(.1,1.1,10)**: Set of resolutions to test the values.
     - **evaluation_metrics=[evaluation.surprise,evaluation.modularity_density]**: Evaluation metrics to be choosen from `CDLIB <https://cdlib.readthedocs.io/en/latest/reference/evaluation.html>`_
     - **kwargs**: Keyword arguments to be sent to `scp.tl.leiden`_.

    **Returns:**

    Adds to adata.obs[key_added] (being key_added="leiden" by default) with the best partition 
    and to adata.uns["leiden"] with the information of the scores fo all resolutions.
    """

    g = adata.uns["neighbors"]["connectivities"]
    G = igraph.Graph(directed=False)
    G.add_vertices(np.unique(g.nonzero()[0]))
    G.add_edges(zip(g.nonzero()[0],g.nonzero()[1]))

    coms = []
    N = {}
    for r in resolutions:
        kwargs["resolution"] = r
        scp.tl.leiden(adata,**kwargs)

        if "key_added" in kwargs.keys():
            key_added=kwargs["key_added"]
        else:
            key_added="leiden"
        
        communities = []
        cluster_ids = adata.obs[key_added].values
        for cluster in np.unique(cluster_ids):
            communities.append(np.where(cluster_ids == cluster)[0])
            
        coms.append(NodeClustering(communities, graph=G, method_name="leiden_"+str(r)))
        N["leiden_"+str(r)] = len(np.unique(cluster_ids))

    rk = evaluation.FitnessRanking(G, coms)
    for metric in evaluation_metrics:
        rk.rank(metric)

    rnk, p_value = rk.friedman_ranking()

    name = np.array([i.alg+"_"+i.param for i in rnk])
    score = np.array([i.score for i in rnk])
    param = np.array([float(i.param) for i in rnk])
    m = {"resolution_best":param[0],"score_best":score[0],"metrics":[str(i) for i in evaluation_metrics]}
    adata.uns["leiden"] = m

    order = np.argsort(param)
    score = score[order]
    param = param[order]
    name = name[order]

    kwargs["resolution"] = adata.uns["leiden"]["resolution_best"]
    scp.tl.leiden(adata,**kwargs)

    adata.uns["leiden"] = m
    adata.uns["leiden"]["scores"] = score
    adata.uns["leiden"]["params"] = param
    adata.uns["leiden"]["n_clusters"] = np.array([N[i] for i in name])
    
    return

def louvain(adata,
    resolutions=np.arange(.1,1.1,.1),
    evaluation_metrics=[evaluation.surprise,evaluation.modularity_density],
    **kwargs):
    """
    Function that constructs over `scp.tl.louvain <https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.louvain.html#scanpy-tl-louvain>`_
    to automatically select the resolution parameter. The optimization algorithm evaluates several metrics for clustering and chooses the best
    resolution among all the metrics using a Friedman Ranking test. The optimization algorithms uses the implementation from 
    `CDLIB <https://cdlib.readthedocs.io/en/latest/index.html>`_.

    **Arguments:**

     - **adata**: Annotated Data object to cluster. It requires "neighbors" to be already computed and present in .uns.
     - **resolutions=np.linspace(.1,1.1,10)**: Set of resolutions to test the values.
     - **evaluation_metrics=[evaluation.surprise,evaluation.modularity_density]**: Evaluation metrics to be choosen from `CDLIB <https://cdlib.readthedocs.io/en/latest/reference/evaluation.html>`_
     - **kwargs**: Keyword arguments to be sent to `scp.tl.louvain`_.

    **Returns:**

    Adds to adata.obs[key_added] (being key_added="louvain" by default) with the best partition 
    and to adata.uns["louvain"] with the information of the scores fo all resolutions.
    """

    g = adata.uns["neighbors"]["connectivities"]
    G = igraph.Graph(directed=False)
    G.add_vertices(np.unique(g.nonzero()[0]))
    G.add_edges(zip(g.nonzero()[0],g.nonzero()[1]))

    coms = []
    N = {}
    for r in resolutions:
        kwargs["resolution"] = r
        scp.tl.louvain(adata,**kwargs)

        if "key_added" in kwargs.keys():
            key_added=kwargs["key_added"]
        else:
            key_added="louvain"
        
        communities = []
        cluster_ids = adata.obs[key_added].values
        for cluster in np.unique(cluster_ids):
            communities.append(np.where(cluster_ids == cluster)[0])
            
        coms.append(NodeClustering(communities, graph=G, method_name="louvain_"+str(r)))
        N["louvain_"+str(r)] = len(np.unique(cluster_ids))

    rk = evaluation.FitnessRanking(G, coms)
    for metric in evaluation_metrics:
        rk.rank(metric)

    rnk, p_value = rk.friedman_ranking()

    name = np.array([i.alg+"_"+i.param for i in rnk])
    score = np.array([i.score for i in rnk])
    param = np.array([float(i.param) for i in rnk])
    m = {"resolution_best":param[0],"score_best":score[0],"metrics":[str(i) for i in evaluation_metrics]}
    adata.uns["louvain"] = m

    order = np.argsort(param)
    score = score[order]
    param = param[order]
    name = name[order]

    kwargs["resolution"] = adata.uns["louvain"]["resolution_best"]
    scp.tl.louvain(adata,**kwargs)

    adata.uns["louvain"] = m
    adata.uns["louvain"]["scores"] = score
    adata.uns["louvain"]["params"] = param
    adata.uns["louvain"]["n_clusters"] = np.array([N[i] for i in name])
    
    return