Workflow Overview
=================

Philosophy
----------

Cellhub is designed to efficiently parallelise the processing of large datasets. Once processed different data slices can be easily extracted directly from the original matrices, aligned and exported for downstream analysis. At the heart of operations is an sqlite database which warehouses the experiment metadata and per-cell statistics.

The workflow can be divided into eight main steps.


1. Quantitation of per-channel libraries
----------------------------------------

The workflow begins with :doc:`pipeline_cellranger_multi.py<pipelines/pipeline_cellranger_multi>`. Input 10X capture channel library identifiers "library_id" and their associated fastq files are specified in the pipeline_cellranger_multi.yml configuration file. The reads from the different libraries will be mapped in parallel.

.. note:: The channel library is considered to be the fundamental "batch" unit of a 10X experiment. Cells captured from the same channel are exposed to the same ambient RNA. Separate genomic libraries are prepared for each 10x channel.

- It is recommend to inspect patterns of ambient RNA using :doc:`pipeline_ambient_rna.py<pipelines/pipeline_ambient_rna>`.

- Per-channel velocity matrices can be prepared using :doc:`pipeline_velocity.py<pipelines/pipeline_velocity>`.

- Cell identification can also be performed with :doc:`pipeline_emptydrops.py<pipelines/pipeline_emptydrops>`.


2. Computation of per-cell statistics
-------------------------------------

Per-cell QC statistics are computed in parallel for each channel library using :doc:`pipeline_cell_qc.py<pipelines/pipeline_cell_qc>`. The pipeline computes various statistics including standard metrics such as percentage of mitochondrial reads, numbers of UMIs and numbers of genes per cell. In addition it can compute scores for custom genesets.

Per-cell celltype predictions are computed in parallel for each channel library using :doc:`pipeline_singleR.py<pipelines/pipeline_singleR>`

.. note:: all per-channel matrices containing computed cell statistics are required to contain "library_id" and "barcode_id" columns.

.. note:: file names of the per-channel matrices are specified as "library_id.tsv.gz" (matrices for different analyses such as e.g. qcmetrics and scrublet scores are written to separate folders).


3. Cell demultiplexing [optional]
---------------------------------

If samples have been multiplexed within channels either genetically or using hash tags a table of barcode_id -> sample_id assignments are prepared using pipeline_demux.py [not yet written].


4. Preparation of the cell database
-----------------------------------

The library and sample metadata, per cell statistics (and demultiplex assignments) are loaded into an sqlite database using :doc:`pipeline_celldb.py<pipelines/pipeline_celldb>`. The pipeline creates a view called "final" which contains the qc and metadata needed for cell selection and downstream analysis.

.. note:: The user is required to supply a tab-separated sample metadata file (e.g. "samples.tsv") via a path in the pipeline_celldb.yml configuration file. It should have columns for library_id, sample_id as well as any other relevant experimental metadata such as condition, genotype, age, replicate, sex etc.


5. Initial assessment of cell quality
-------------------------------------

This is performed manally by the data analyst as it is highly dataset-specific. Per cell QC statistics can be easily retrieved from the celldb for plotting along with singleR scores from the cellhub API.


6. Fetching of cells for downstream analysis
--------------------------------------------

Cells are fetched using :doc:`pipeline_fetch_cells.py<pipelines/pipeline_fetch_cells>`. The user specifies the cells that they wish to retrieve from the "final" table (see step 4) via an sql statement in the pipeline_fetch_cells.yml configuration file. The pipeline will extract the cells and metadata from the original matrices and combine them into market matrices and anndata objects for downstream analyses.

It is recommended to fetch cells into a new directory. By design fetching of a single dataset per-directory is supported.

The pipeline supports fetching of velocity information.

.. note:: The retrieved metadata will include a "sample_id" column. From this point onwards it is natural to think of the "sample_id" as the unit of interest. The "library_ids" remain in the metadata along with all the qc statistics to facilitate downstream investigation of batch effects and cell quality.

7. ADT normalization [optional]
-------------------------------
If samples included the ADT modality, :doc:`pipeline_adt_norm.py<pipelines/pipeline_adt_norm>` normalizes the antibody counts for the high-quality fetched cells in the previous step. Normalized ADT can be then used for downstream integration. The pipeline implements 3 normalization methodologies: DSB, median-based, and CLR. The user can specify the feature space.

8. Integration
--------------

Integration of samples is performed manually by the user because it is highly dataset specific. Different integration algorithms are needed for different contexts. Strategies for HVG selection and modelling of covariates need to be considered by the data analyst on a case by case basis.

9. Clustering analysis
----------------------

Clustering analysis is performed with  :doc:`pipeline_cluster<pipelines/pipeline_cluster>`.


Workflow Diagram
================

The diagram is now a little out of date with respect to configuration of the pipeline inputs but provides a useful depiction of the overall workflow.

.. image:: images/cellhub-devel-schema.png
