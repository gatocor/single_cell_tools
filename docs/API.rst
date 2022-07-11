#################
API
#################

===========================
Tools
===========================

.. automodule:: single_cell_tools.tools
    :members: qc_metrics, filter_cells, rank_genes_to_excel, enrichr_to_excel, downsample, common_genes, scmap_annotate, scmap_projection

===========================
Hyperparameter selection
===========================

---------------------------
PCA
---------------------------
.. automodule:: single_cell_tools.hyperparameter_selection
    :members: pca_explained_variance, pca_permutation_test

---------------------------
Clustering
---------------------------
.. automodule:: single_cell_tools.hyperparameter_selection
    :members: leiden, louvain

===========================
Evaluation
===========================
.. automodule:: single_cell_tools.evaluation
    :members: kBET

===========================
Plotting
===========================
.. automodule:: single_cell_tools.plot
    :members: plot_base, vline, hline, plot_scatter_labels, plot_scatter_lines, plot_scatter_pies
