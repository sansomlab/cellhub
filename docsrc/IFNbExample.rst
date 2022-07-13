IFNb PBMC example
=================

Setting up
----------

1. Clone the example template to a local folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the folders and files for the example into a local folder.

  cp -r /path/to/cellhub/examples/ifnb_pbmc/* .

This will create 2 folders:

- "cellhub" where the preprocessing pipelines will be run and where the cellhub database will be created
- "integration" where integration is to be performed, an example integrated anndata object is provided.
- "cluster" where we can perform downstream analysis with "cellhub cluster".


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

  cellhub cell_qc -v5 -p20 make full


4. Running emptydrops and investigating ambient RNA (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If desired we can run emptydrops::

  # cellhub emptydrops config

  cellhub emptydrops -v5 -p20 make full

And investigate the ambient rna::

  # cellhub ambient_rna config

  cellhub ambient_rna -v5 -p20 make full


5. Loading the cell statistics into the celldb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cell QC statistics and metadata ("samples.tsv") are next loaded into a local sqlite database::

  # cellhub celldb config

  cellhub celldb -v5 -p20 make full


6. Run pipeline_singleR
^^^^^^^^^^^^^^^^^^^^^^^^

Single R is run on all the cells so that the results are avaliable to help with QC
as well as downstream analysis.

  # cellhub singleR config
  
  cellhub singleR -v5 -p20 make full.
  
As noted: :doc:`in the pipeline_singleR inputs section <pipelines/pipeline_singleR>` the celldex references
neede to be stashed before the pipeline is run.


7. Run pipeline_annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^

This short pipeline retrieves Ensembl and KEGG annotations needed for downstream analysis.

  # cellhub annotation config
  
  cellhub annotation -v5 -p10 make full
  
Please note that the specified Ensembl version should match that used for the cellranger reference trancriptome.


Performing cell QC
------------------


8. Assessment of cell quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is left to the reader to perform manually because it needs to be carefully tailored to individual datasets.


Performing downstream analysis
------------------------------


9. Fetch cells for integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use pipeline_fetch_cells to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are 
specified in the pipeline_fetch_cells.yml) file::

  # It is recommended to fetch the cells in to a seperate directory for integration.
  cd ../integration

  # cellhub fetch_cells config
  cellhub fetch_cells  -v5 -p20 make full


10. Integration
^^^^^^^^^^^^^^

This step is performed manually because it is highly dataset specific. Different integration algorithms are needed for different contexts and strategies for HVG selection and modelling of covariates need
to be considered by the data analyst on a case by case basis.

The result should be saved as an anndata file as described: :doc:`in the pipeline_cluster inputs section <pipelines/pipeline_cluster>`.


11. Clustering analysis
^^^^^^^^^^^^^^^^^^^^^^^

Cluster analysis is performed with pipeline cluster. It is recommedend to do this in a seperate directory.

  # change into a new directory
  cd ../cluster.dir

  # checkout a the yml file and configure the pipeline with the location of the cellhub directory and
  # integrated anndata object
  cellhub cluster config
  
  # a suitable yml file has been provided so we can now launch the pipeline
  cellhub cluster -v5 -p200 make full
