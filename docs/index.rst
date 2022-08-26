.. cellhub documentation master file, created by
   sphinx-quickstart on Tue May  4 15:09:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cellhub: a single-cell data analysis superhighway
=================================================

Cellhub provides an end-to-end scaleable workflow for the pre-processing, warehousing and analysis of data from millions of single-cells. It aims to brings together best practice solutions, including for read alignment (`Cellranger <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger`_>), ambient RNA correction (`CellBender <https://github.com/broadinstitute/CellBender>`_), de-multiplexing and de-hashing (`GMMDemux <https://github.com/CHPGenetics/GMM-Demux>`_), cell type prediction (`singleR <https://bioconductor.org/packages/release/bioc/html/SingleR.html>`_) and cluster analysis (`Scanpy <https://scanpy.readthedocs.io/en/stable/>`_) into a cohesive set of easy to use analysis pipelines. It relies on the `cgat-core workflow management system <https://github.com/cgat-developers/cgat-core>`_ to leverage the power of high-performance compute clusters for high-throughput parallel processing. Pipelines cross-talk via a defined API allowing for easily extension or modification of the workflow. Cells and their associated qc-statistics and metadata are indexed in a central SQLite database from which arbitrary subsets can be easily fetched for downstream analysis in `anndata format <https://anndata.readthedocs.io/en/latest/>`_. The clustering pipeline allows rapid evaluation of different pre-processing and integration strategies at different clustering resolutions, through generation of pdf reports and `cellxgene <https://github.com/chanzuckerberg/cellxgene>` objects.

Cellhub was originally developed to support the single-cell component of the `University of Oxford's COMBAT COVID-19 project<https://doi.org/10.1016/j.cell.2022.01.012>`_. 

Cellhub is currently alpha software. More detailed examples and tutorials will follow soon.

.. toctree::
   :maxdepth: 2

   overview.rst
   installation.rst
   usage.rst
   examples.rst
   pipelines.rst
   api.rst
   tasks.rst
   contributing.rst

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
