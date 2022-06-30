def qc_metrics(adata, mtgenes):
    """
    Function that computes the basic QC metrics (#counts per cell, #expressed_genes per cell, mitochondrial fraction)

    **Arguments:**
     - **adata**: Annotated data where to make the metrics.
     - **mtgenes**: Bool array of the length of .var with True in the genes that are mitochondrial.

    **Returns:**
     Annotated data with #counts, #expressed_genes, mtFraction added to .obs.
    """
    adata.obs["#counts"] = adata.X.sum(axis=1)
    adata.obs["#expressed_genes"] = (adata.X > 0).sum(axis=1)
    adata.obs["mtFraction"] = (adata[:,mtgenes].X).sum(axis=1)
    adata.obs["mtFraction"] /= adata.obs["#counts"].values

    return