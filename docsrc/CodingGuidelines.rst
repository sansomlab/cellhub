Coding Guidelines
=================

Repository layout
-----------------

.. list-table:: Repository layout
   :widths: 25 100
   :header-rows: 1

   * - Folder
     - Contents
   * - cellhub
     - The cellhub Python module which contains the set of CGAT-core pipelines
   * - cellhub/tasks
     - The cellhub tasks submodule which contains helper functions and pipeline task definitions
   * - cellhub/reports
     - The latex source files used for building the reports
   * - Python
     - Python worker scripts that are executed by pipeline tasks
   * - R/cellhub
     - The R cellhub library
   * - R/scripts
     - R worker scripts that are executed by pipeline tasks
   * - docsrc
     - The documentation source files in restructured text format for compilation with sphinx
   * - docs
     - The compiled documentation
   * - examples
     - Configuration files for example datasets
   * - conda
     - Conda environment, TBD
   * - tests
     - TBD


Coding style
------------

Currently we are working to improve and standardise the coding style:

Python
^^^^^^

* Python code should be `pep8 <https://www.python.org/dev/peps/pep-0008/>`_ compliant. Compliance checks will be enforced soon.

* Arguments to Python scripts should be parsed with argparse.

* Logging in Python scripts should be performed with the standard library "logging" module, written to stdout and redirected to a log file in the pipeline task.

R
^

* R code should follow the `tidyverse style guide <https://style.tidyverse.org>`_. Please do not use right-hand assignment.

* Arguments to R scripts should be parsed with optparse.

* Logging in R scripts should be performed with the "futile.logger" library, written to stdout and redirect to a log file specified in the pipeline task.

* Otherwise, to write to stdout from R scripts use message() or warning(). Do not use print() or cat().


Writing pipelines
-----------------

The pipelines live in the "cellhub" python module.

Auxiliary task functions live in the "cellhub/task" python sub-module.

.. note:: Tasks of more than a few lines should be abstracted into appropriately named sub-modules.

In the notes below "xxx" denotes the name of a pipeline such as e.g. "cell_qc".

1. Paths should never be hardcoded in the pipelines - rather they must be read from the yaml files.
2. Yaml configuration files should be named pipeline_xxx.yml
3. The output of individual pipelines should be written to a subfolder name "xxx.dir" to keep the root directory clean (it should only contain these directories and the yml configuration files!).
4. Pipelines that generate cell-level information for down-stream analysis must read their inputs from the api and register their public outputs to the API, see :doc:`API<API>`. If you need information from an upstream pipeline that is not present on the API please raise an issue.
5. We are documenting the pipelines using the sphinx "autodocs" module so please maintain informative rst docstrings.


Cell indexing
-------------

* Tables of per-cell information registered on the API must have columns "barcode" (for the original Cellranger assigned barcode, "-1" suffix not removed) and "library_id". These two columns are used by pipeline_celldb.py to join the tables in the database.

* Downstream of fetch_cells, when cells are aggregated across libraries, we use unique cell identifiers made by combining the barcode and library id in [AGCT]-library_id format (with the "-1" suffix now removed from the original barcode). The unique cell identifiers are used to populate the anndata.obs.index and the "barcode_id" column in tsv files where needed.



Yaml configuration file naming
------------------------------

The cgat-core system only supports configuration files name "pipeline.yml".

We work around this by overriding the cgat-core functionality using a helper function in cellhub.tasks.control as follows::

  import Pipeline as P
  import cellhub.tasks.control as C

  # Override function to collect config files
  P.control.write_config_files = C.write_config_files

Default yml files must be located at the path cellhub/yaml/pipeline_xxx.yml


Writing and compiling the documentation
---------------------------------------

The source files for the documentation are found in::

  docsrc

The documentation source files for the pipelines can be found in::

  docsrc/pipelines

To build the documentation cd to the docsrc folder and run::

  make github

This will build the documentation and copy the latex output to the "docs" folder. You then need to cd to the "docs" folder and run::

  make

To compile the pdf.

When the repo is made public we will switch to using html documentation on readthedocs. Unfortunately there is no straightforward solution for private html hosting.
