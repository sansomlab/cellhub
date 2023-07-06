10k Human PBMCs 5' v2.0 example
===============================

The 10x `10k Human PBMCs, 5' v2.0 <https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-5-v-2-0-chromium-controller-2-standard-6-1-0>`_ dataset contains multimodal GEX, TCR and BCR data from a single human subject.

Example configuration files ("samples.tsv", "libraries.tsv" and "pipeline_cellranger.yml") for processing this dataset are provided in the "examples/10k_pbmc" folder.

.. Note::  Due to a quirk with "cellranger vdj" fastq data for vdj libraries from different flow cells (for the same sample) needs to be presented in different folders. In other words, if fastqs for a vdj sample from different flow cells are presented in the same folder with different "sample" prefix the "cellranger vdj" command fails.

For downstream analysis please see the :doc:`IFNb example</examples/ifnb>`.

