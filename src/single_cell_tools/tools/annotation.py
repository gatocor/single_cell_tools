import pandas as pd
import numpy as np
import json
import requests

def rank_genes_to_excel(adata,gene_name,filename,n_genes,layer="Raw"):
    """
    Handy function to write the rank_genes_groups from scanpy to an excel table.

    **Arguments:**\n
     - **adata**: Annotated Data with the scanpy.tl.rank_genes_groups runned.
     - **gene_name**: Columns in .var to use for recognising the gene.
     - **filename**: Name of file where to save the results with ".xlsx" in the end.
     - **n_genes**: List of the first ranked genes of each cluster to put in the list. 

    **Returns:**

    Excel file with the DE genes in different pages of the file.
    """
            
    groupby = adata.uns["rank_genes_groups"]["params"]["groupby"]
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for k,group in enumerate(adata.uns["rank_genes_groups"]["names"].dtype.names):
                
        l = pd.DataFrame(columns=["names","logfoldchanges","pvals","pvals_adj","scores"])
        l.loc[:,"names"] = adata.var.loc[[n[k] for n in adata.uns["rank_genes_groups"]["names"]][:n_genes],gene_name].values
        pop1 = adata.obs.loc[:,groupby]==group
        pop2 = adata.obs.loc[:,groupby]!=group
        l.loc[:,"logfoldchanges"] = [np.log(adata[pop1,adata.var.loc[:,gene_name]==j].layers[layer].mean()/adata[pop2,adata.var.loc[:,gene_name]==j].layers[layer].mean()) for j in l.loc[:,"names"]]
        l.loc[:,"mean_expression_in_group"] = [adata[pop1,adata.var.loc[:,gene_name]==j].layers[layer].mean() for j in l.loc[:,"names"]]
        l.loc[:,"pvals"] = [n[k] for n in adata.uns["rank_genes_groups"]["pvals"]][:n_genes]
        l.loc[:,"pvals_adj"] = [n[k] for n in adata.uns["rank_genes_groups"]["pvals_adj"]][:n_genes]
        l.loc[:,"scores"] = [n[k] for n in adata.uns["rank_genes_groups"]["scores"]][:n_genes]
        
        l.to_excel(writer, sheet_name="cluster_"+str(group))
            
    writer.save()

    return

def enrichr_to_excel(adata,library,filename,n_genes=200):
    """
    Function that sends to `Enrichr <https://maayanlab.cloud/Enrichr/>`_ the genes computed by scp.tl.rank_genes_group, 
    retrieves the information from one of the libraries and saves the information in an excel file with each cluster as excel.

    **Arguments:**\n
     - **adata**: Annotated Data with the scanpy.tl.rank_genes_groups runned.
     - **library**: Library from `enrich libraries <https://maayanlab.cloud/Enrichr/#libraries>`_. 
     - **filename**: Name of file where to save the results with ".xlsx" in the end.
     - **n_genes**: List of the first ranked genes of each cluster to put in the list. 

    **Returns:**

    Excel file with the Enriche queries in different pages of the file.
    """

    ENRICHR_URL_ADD = 'https://maayanlab.cloud/Enrichr/addList'
    description = 'Example gene list'

    ENRICHR_URL_QUERY = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for i,group in enumerate(adata.uns["rank_genes_groups"]["names"].dtype.names):
        #Send the list
        genes = adata.var.loc[[n[i] for n in adata.uns["rank_genes_groups"]["names"]][:n_genes],"gene_name"].values
        genes_str = '\n'.join(genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        user_list_id = json.loads(response.text)["userListId"]

        #Get the results of the library
        gene_set_library = library
        response = requests.get(
            ENRICHR_URL_QUERY + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)[gene_set_library]

        data = pd.DataFrame(data)
        data.columns = ["Index","Name","P-value","OddsRatio","Combined Score","Genes","Adjusted P-value","Old P-value","Old Adjusted P-value"]
        data.to_excel(writer,sheet_name="Cluster_"+str(group))

    writer.save()  

def enrichr_to_excel_old(adata,library,filename,n_genes=200):
    """
    Function that sends to `Enrichr <https://maayanlab.cloud/Enrichr/>`_ the genes computed by scp.tl.rank_genes_group, 
    retrieves the information from one of the libraries and saves the information in an excel file with each cluster as excel.

    **Arguments:**\n
     - **adata**: Annotated Data with the scanpy.tl.rank_genes_groups runned.
     - **library**: Library from `enrich libraries <https://maayanlab.cloud/Enrichr/#libraries>`_. 
     - **filename**: Name of file where to save the results with ".xlsx" in the end.
     - **n_genes**: List of the first ranked genes of each cluster to put in the list. 

    **Returns:**

    Excel file with the Enriche queries in different pages of the file.
    """

    ENRICHR_URL_ADD = 'https://maayanlab.cloud/Enrichr/addList'
    description = 'Example gene list'

    ENRICHR_URL_QUERY = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for i,group in enumerate(adata.obs[adata.uns["rank_genes_groups"]["params"]["groupby"]].unique()):
        #Send the list
        genes = adata.var.loc[[n[i] for n in adata.uns["rank_genes_groups"]["names"]][:n_genes],"gene_name"].values
        genes_str = '\n'.join(genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(ENRICHR_URL_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        user_list_id = json.loads(response.text)["userListId"]

        #Get the results of the library
        gene_set_library = library
        response = requests.get(
            ENRICHR_URL_QUERY + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)[gene_set_library]

        data = pd.DataFrame(data)
        data.columns = ["Index","Name","P-value","OddsRatio","Combined Score","Genes","Adjusted P-value","Old P-value","Old Adjusted P-value"]
        data.to_excel(writer,sheet_name="Cluster_"+str(group))

    writer.save()  
