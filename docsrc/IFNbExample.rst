IFNb PBMC example
=================

To get started create a suitable directory and cd into it.

Copy the sample_metadata.tsv file from the examples/ifnb_pbmc/sample_metadata.tsv [TODO!]


1. Running cellranger
=====================

The first step is to configure and run pipeline_cellranger_multi::

  python /path/to/cellhub-devel/pipelines/pipeline_cellranger_multi.py config

  python /path/to/cellhub-devel/pipelines/pipeline_cellranger_multi.py make full -v5 -p20


2. Running the cell qc pipeline
===============================

Next we run the cell qc pipeline::

  python /path/to/cellhub-devel/pipelines/pipeline_cell_qc.py config

  python /path/to/cellhub-devel/pipelines/pipeline_cell_qc.py make full -v5 -p20


3. Running emptydrops (optional)
================================

If desired we can run emptydrops::

  python /path/to/cellhub-devel/pipelines/pipeline_emptydrops.py config

  python /path/to/cellhub-devel/pipelines/pipeline_emptydrops.py make full -v5 -p20


4. Loading the cell statistics into the celldb
==============================================

The cell QC statistics and metadata are next loaded into a local sqlite database::

  python /path/to/cellhub-devel/pipelines/pipeline_celldb.py config

  python /path/to/cellhub-devel/pipelines/pipeline_celldb.py make full -v5 -p20


5. Fetching cells for downstream analysis
=========================================

We use pipeline_fetch_cells.py to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are specified in the pipeline_fetch_cells.yml) file::

  python /path/to/cellhub-devel/pipelines/pipeline_fetch_cells.py config

  python /path/to/cellhub-devel/pipelines/pipeline_fetch_cells.py make full -v5 -p20


6. Integration
==============

pipeline_integration.py supports integration of the data with harmony, bbknn and scanorama::

  python /path/to/cellhub-devel/pipelines/pipeline_integration.py config

  python /path/to/cellhub-devel/pipelines/pipeline_integration.py make full -v5 -p20


7. Export for downstream analysis
=================================

Finally we can export the integrated anndata object to e.g. a Seurat object for downstream analysis::

  python /path/to/cellhub-devel/pipelines/pipeline_export.py config

  python /path/to/cellhub-devel/pipelines/pipeline_export.py make full -v5 -p20
