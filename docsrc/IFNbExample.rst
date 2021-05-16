IFNb PBMC example
=================

Setting up
----------

1. Clone the example template to a local folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the folders and files for the example into a local folder.

  cp -r /path/to/cellhub/examples/ifnb_pbmc/* .

This will create 3 folders:

- "cellhub" where the preprocessing pipelines will be run and where the cellhub database will be created
- "integration" where data for cells of interest (e.g. passing qc) will be fetched and integrated
- "scxl" where we can perform downstream analysis with `pipeline_scxl.py <https://github.com/sansomlab/tenx>`_.


Running the pre-processing pipelines and creating the database
--------------------------------------------------------------

2. Running Cellranger
^^^^^^^^^^^^^^^^^^^^^

The first step is to configure and run pipeline_cellranger_multi. We already have a pre-configured yml file so we can skip this step but the syntax is included for reference here and also for the other steps: ::

  # enter the cellhub directory

  cd cellhub

  # cellhub cellranger_multi config

Edit the pipeline_cellranger_multi.yml file as appropriate to point to folders containing fastq files extracted from the original BAM files submitted by `Kang et. al. <https://doi.org/10.1038/nbt.4042>`_ to GEO (GSE96583). The GEO identifiers are: unstimulated (GSM2560248) and stimulated (GSM2560249). The fastqs can be extracted with the `10X bamtofastq tool <https://support.10xgenomics.com/docs/bamtofastq>`_.

We run the pipeline as follows: ::

  cellhub cellranger_multi make full -v5 -p20

If you have not run "python setup.py devel" pipelines can instead be launched directly. In this case the equivalent command would be::

  python path/to/cellhub-devel/cellhub/pipeline_cellranger_multi.py make full -v5 -p20 --pipeline-log=pipeline_cellranger_multi.py

.. note:: when launching pipelines directly if the "--pipeline-log" parameter is not specified the log file will be written to "pipeline.log".


3. Running the cell qc pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next we run the cell qc pipeline::

  # cellhub cell_qc config

  cellhub cell_qc make full -v5 -p20


4. Running emptydrops and investigating ambient RNA (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If desired we can run emptydrops::

  # cellhub emptydrops config

  cellhub emptydrops make full -v5 -p20

And investigate the ambient rna::

  # cellhub ambient_rna config

  cellhub ambient_rna make full -v5 -p20


5. Loading the cell statistics into the celldb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cell QC statistics and metadata ("samples.tsv") are next loaded into a local sqlite database::

  # cellhub celldb config

  cellhub celldb make full -v5 -p20


Performing downstream analysis
------------------------------


6. Fetching cells for downstream analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


We use pipeline_fetch_cells to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are specified in the pipeline_fetch_cells.yml) file::

  # change into the integration directory
  cd ../integration

  # cellhub fetch_cells config
  cellhub fetch_cells make full -v5 -p20



7. Integration
^^^^^^^^^^^^^^

Pipeline_integration supports integration of the data with harmony, bbknn and scanorama. In this example the location of the data is specified in the "integration.tsv" file as per the path given in the "pipeline_integration.yml" file. ::

  # cellhub integration config

  cellhub integration make full -v5 -p20

.. warning:: pipeline_integration.py will be moving to a new sansomlab/scxl repository (along with pipeline_scxl from sansomlab/tenx).


8. Export for downstream analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally we can export the integrated anndata object to e.g. a Seurat object for downstream analysis::

  # cellhub export config

  cellhub export make full -v5 -p20

.. warning:: pipeline_export.py will be moving to the new sansomlab/scxl repository


9. Clustering analysis with pipeline_scxl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cluster analysis can be performed with `pipeline_scxl.py <https://github.com/sansomlab/tenx>`_. A suitable configuration file for working with the harmony aligned seurat object is provided in the examples/infb_pbmc/scxl/ folder::

  # change into the sxcl directory
  cd ../scxl.dir

  # link in the exported Seurat object from step 8
  mkdir integrated.seurat.dir
  ln -s ../integration/export.dir/pbmc.exp.dir/seurat_object.rds integrated.seurat.dir/begin.rds

  # a suitable yml file has been provided so we can now launch the pipeline
  python /path/to/tenx/pipelines/pipeline_scxl.py make full -v5 -p200
