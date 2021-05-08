IFNb PBMC example
=================

To get started create a suitable directory and cd into it.

Copy the yml, sample.tsv, integration.tsv and export.tsv configuration and metadata files for this example to the working folder::

  cp /path/to/cellhub/examples/ifnb_pbmc/* .


1. Running cellranger
=====================

The first step is to configure and run pipeline_cellranger_multi. We already have a pre-configured yml file so we can skip this step but the syntax is included for referencehere and also for the other steps ::

  # python /path/to/cellhub-devel/pipelines/pipeline_cellranger_multi.py config

Edit the pipeline_cellranger_multi.yml file as appropriate to point to folders containing fastq samples extracted from the original BAM files submitted by `Kang et. al. <https://doi.org/10.1038/nbt.4042>` to GEO (GSE96583). The GEO identifiers are: unstimulated (GSM2560248) and stimulated (GSM2560249). The fastqs can be extracted with the `10X bamtofastq tool <https://support.10xgenomics.com/docs/bamtofastq>`.

And run the pipeline being careful to redirect the log to an appropriately named file.::

  python /path/to/cellhub-devel/pipelines/pipeline_cellranger_multi.py make full -v5 -p20 --pipeline-logfile=pipeline_cellranger_multi.log


2. Running the cell qc pipeline
===============================

Next we run the cell qc pipeline::

  # python /path/to/cellhub-devel/pipelines/pipeline_cell_qc.py config

  python /path/to/cellhub-devel/pipelines/pipeline_cell_qc.py make full -v5 -p20 --pipeline-logfile=pipeline_cell_qc.log


3. Running emptydrops and investigating ambient RNA (optional)
==============================================================

If desired we can run emptydrops::

  # python /path/to/cellhub-devel/pipelines/pipeline_emptydrops.py config

  python /path/to/cellhub-devel/pipelines/pipeline_emptydrops.py make full -v5 -p20 --pipeline-logfile=pipeline_emptydrops.log

And investigate the ambient rna::

  # python /path/to/cellhub-devel/pipelines/pipeline_ambient_rna.py config

  python /path/to/cellhub-devel/pipelines/pipeline_ambient_rna.py make full -v5 -p20 --pipeline-logfile=pipeline_emptydrops.log


4. Loading the cell statistics into the celldb
==============================================

The cell QC statistics and metadata are next loaded into a local sqlite database::

  # python /path/to/cellhub-devel/pipelines/pipeline_celldb.py config

  python /path/to/cellhub-devel/pipelines/pipeline_celldb.py make full -v5 -p20 --pipeline-logfile=pipeline_celldb.log


5. Fetching cells for downstream analysis
=========================================

We use pipeline_fetch_cells.py to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are specified in the pipeline_fetch_cells.yml) file::

  # python /path/to/cellhub-devel/pipelines/pipeline_fetch_cells.py config

  python /path/to/cellhub-devel/pipelines/pipeline_fetch_cells.py make full -v5 -p20 --pipeline-logfile=pipeline_fetch_cells.log


6. Integration
==============

pipeline_integration.py supports integration of the data with harmony, bbknn and scanorama::

  # python /path/to/cellhub-devel/pipelines/pipeline_integration.py config

  python /path/to/cellhub-devel/pipelines/pipeline_integration.py make full -v5 -p20 --pipeline-logfile=pipeline_integration.log


7. Export for downstream analysis
=================================

Finally we can export the integrated anndata object to e.g. a Seurat object for downstream analysis::

  # python /path/to/cellhub-devel/pipelines/pipeline_export.py config

  python /path/to/cellhub-devel/pipelines/pipeline_export.py make full -v5 -p20 --pipeline-logfile=pipeline_export.log


8. Perform downstream analysis
==============================

Downstream analysis can be peformed with `pipeline_scxl.py <https://github.com/sansomlab/tenx>`. A suitable configuration file for working with the harmony aligned seurat object is provided in the examples/infb_pbmc/scxl/ folder::

  mkdir scxl.dir
  cd scxl.dir
  mkdir integrated.seurat.dir
  ln -s ../export.dir/pbmc.exp.dir/seurat_object.rds integrated.seurat.dir/begin.rds
  cp /path/to/cellhub/examples/ifnb_pbmc/pipeline_scxl/pipeline.yml .
  python /path/to/tenx/pipelines/pipeline_scxl.py make full -v5 -p200
