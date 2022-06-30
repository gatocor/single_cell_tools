#################
API
#################

=================
Pipeline
=================

.. automodule:: single_cell_tools.pipeline
    :members:

===========================
Tools
===========================

.. automodule:: single_cell_tools.tools
    :members: qc_metrics, filter_cells

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
    :members: batch_kBET

===========================
Plotting
===========================
.. automodule:: single_cell_tools.plot
    :members: plot_base, vline, hline
