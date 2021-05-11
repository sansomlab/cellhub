IFNb PBMC example
=================

To get started create a suitable directory and cd into it.

Copy the yml, sample.tsv, integration.tsv and export.tsv configuration and metadata files for this example to the working folder::

  cp /path/to/cellhub/examples/ifnb_pbmc/* .


1. Running Cellranger
---------------------

The first step is to configure and run pipeline_cellranger_multi. We already have a pre-configured yml file so we can skip this step but the syntax is included for reference here and also for the other steps: ::

  # cellhub cellranger_multi config

Edit the pipeline_cellranger_multi.yml file as appropriate to point to folders containing fastq files extracted from the original BAM files submitted by `Kang et. al. <https://doi.org/10.1038/nbt.4042>`_ to GEO (GSE96583). The GEO identifiers are: unstimulated (GSM2560248) and stimulated (GSM2560249). The fastqs can be extracted with the `10X bamtofastq tool <https://support.10xgenomics.com/docs/bamtofastq>`_.

We run the pipeline as follows: ::

  cellhub cellranger_multi make full -v5 -p20

If you have not run "python setup.py devel" pipelines can instead be launched directly. In this case the equivalent command would be::

  python path/to/cellhub-devel/cellhub/pipeline_cellranger_multi.py make full -v5 -p20 --pipeline-log=pipeline_cellranger_multi.py

.. note:: when launching pipelines directly if the "--pipeline-log" parameter is not specified the log file will be written to "pipeline.log".


2. Running the cell qc pipeline
-------------------------------

Next we run the cell qc pipeline::

  # cellhub cell_qc config

  cellhub cell_qc make full -v5 -p20


3. Running emptydrops and investigating ambient RNA (optional)
--------------------------------------------------------------

If desired we can run emptydrops::

  # cellhub emptydrops config

  cellhub emptydrops make full -v5 -p20

And investigate the ambient rna::

  # cellhub ambient_rna config

  cellhub ambient_rna make full -v5 -p20


4. Loading the cell statistics into the celldb
----------------------------------------------

The cell QC statistics and metadata ("samples.tsv") are next loaded into a local sqlite database::

  # cellhub celldb config

  cellhub celldb make full -v5 -p20


5. Fetching cells for downstream analysis
-----------------------------------------

We use pipeline_fetch_cells to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are specified in the pipeline_fetch_cells.yml) file::

  # cellhub fetch_cells config

  cellhub fetch_cells make full -v5 -p20


6. Integration
--------------

Pipeline_integration supports integration of the data with harmony, bbknn and scanorama. In this example the location of the data is specified in the "integration.tsv" file as per the path given in the "pipeline_integration.yml" file. ::

  # cellhub integration config

  cellhub integration make full -v5 -p20


7. Export for downstream analysis
---------------------------------

Finally we can export the integrated anndata object to e.g. a Seurat object for downstream analysis::

  # cellhub export config

  cellhub export make full -v5 -p20


8. Perform downstream analysis
------------------------------

Downstream analysis can be performed with `pipeline_scxl.py <https://github.com/sansomlab/tenx>`_. A suitable configuration file for working with the harmony aligned seurat object is provided in the examples/infb_pbmc/scxl/ folder::

  mkdir scxl.dir
  cd scxl.dir
  mkdir integrated.seurat.dir
  ln -s ../export.dir/pbmc.exp.dir/seurat_object.rds integrated.seurat.dir/begin.rds
  cp /path/to/cellhub/examples/ifnb_pbmc/pipeline_scxl/pipeline.yml .
  python /path/to/tenx/pipelines/pipeline_scxl.py make full -v5 -p200
