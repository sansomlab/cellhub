IFNb PBMC example
=================

Setting up
----------

1. Clone the example template to a local folder
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the folders and files for the example into a local folder: ::

  cp -r /path/to/cellhub/examples/ifnb_pbmc/* .

This will create 3 folders:

- "cellhub" where the preprocessing pipelines will be run and where the cellhub database will be created
- "integration" where integration is to be performed.
- "cluster" where we can perform downstream analysis with "cellhub cluster".

The folders contain the necessary configuration files: please edit the cellhub/libraries.tsv file to point to the location of the FASTQ files on your system.

Running the pre-processing pipelines and creating the database
--------------------------------------------------------------

2. Running Cellranger
^^^^^^^^^^^^^^^^^^^^^

The first step is to configure and run pipeline_cellranger. This pipeline takes three inputs (i) Information about the biological samples, number of cells expected and the 10x chemistry version are specified in a "samples.tsv" file. (ii) Input 10X  channel library identifiers "library_id", sample prefixes, library types and FASTQ paths are specified via a tab delimited "libraries.tsv" file. (iii) A pipeline_cellranger.yml file is used to configure general options such as computional resource specification and the location of the genomic references. For more details please see: :doc:`pipeline_cellranger.py</pipelines/pipeline_cellranger>`. 

For this example, preconfigured "samples.tsv", "libraries.tsv" and "pipeline_cellranger.yml" files are provided in the cellhub directory. 

Edit the "libraries.tsv" file "fastq_path" column as appropriate to point to folders containing fastq files extracted from the original BAM files submitted by `Kang et. al. <https://doi.org/10.1038/nbt.4042>`_ to GEO (GSE96583). The GEO identifiers are: unstimulated (GSM2560248) and stimulated (GSM2560249). The fastqs can be extracted with the `10X bamtofastq tool <https://support.10xgenomics.com/docs/bamtofastq>`_.

The pipeline is run as follows: ::

  cellhub cellranger make full -v5 -p20

Finally, the count matrices must be manually registered on the API for downstream analysis: ::

  cellhub cellranger make useCounts

.. note:: When processing other datasets the "samples.tsv" and "libraries.tsv" files must be created by the user. For more details on constructing these files please see :doc:`pipeline_cellranger.py</pipelines/pipeline_cellranger>`. A template 'pipeline_cellranger.yml' file can be obtained using the "config" command which is common to all cellhub piplines.:

  # cellhub cellranger config


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


5. Run pipeline_singleR
^^^^^^^^^^^^^^^^^^^^^^^^

Single R is run on all the cells so that the results are avaliable to help with QC
as well as downstream analysis::

  # cellhub singleR config
  
  cellhub singleR make full -v5 -p20
  
As noted: :doc:`in the pipeline_singleR inputs section <pipelines/pipeline_singleR>` the celldex references
need to be stashed before the pipeline is run.


6. Loading the cell statistics into the celldb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cell QC statistics and metadata ("samples.tsv") are next loaded into a local sqlite database::

  # cellhub celldb config

  cellhub celldb make full -v5 -p20


7. Run pipeline_annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^

This  pipeline retrieves Ensembl and KEGG annotations needed for downstream analysis.::

  # cellhub annotation config
  
  cellhub annotation make full -v5 -p10 
  
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

We use pipeline_fetch_cells to retrieve the cells we want for downstream analysis. (QC thresholds and e.g. desired samples are specified in the pipeline_fetch_cells.yml) file::

  # It is recommended to fetch the cells in to a seperate directory for integration.
  cd ../integration

  # cellhub fetch_cells config
  cellhub fetch_cells make full -v5 -p20


10. Integration
^^^^^^^^^^^^^^^

Run the provided jupyter notebook to perform a basic Harmony integration of the data and to save it in the appropriate anndata format (see :doc:`in the pipeline_cluster inputs section <pipelines/pipeline_cluster>`) is provided.


11. Clustering analysis
^^^^^^^^^^^^^^^^^^^^^^^

Cluster analysis is performed with pipeline cluster (a seperate directory is recommended for this so that multiple clustering runs can be performed as required).::

# change into the clustering directory
cd ../cluster

# A suitable `pipeline_cluster.yml` has already been provided in the `cluster` directory.

# If needed, you can generate a new one with:
# cellhub cluster config
# (Note: This will default to `rdim_name: scVI`, which may not match Harmony output)

cellhub cluster make full -v5 -p200


The pdf reports and excel files generated by the pipeline can be found in the "reports.dir" subfolder.

For interactive visulation, the results are provided in cellxgene format. To view the cellxgene.h5ad files, you will first need toinstall cellxgene with "pip install cellxgene". The cellxgene viewer can then be launched with: ::

  # substitute "{x}" with the number integrated components used for the clustering run.
  cellxgene --no-upgrade-check launch out.{x}.comps.dir/cellxgene.h5ad
