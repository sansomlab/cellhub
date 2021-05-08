Usage
=====


Configuring and running pipelines
---------------------------------

Run the cellhub --help command to view the help documentation and find available pipelines to run cellhub.

The cellhub pipelines are written using `cgat-core <https://github.com/cgat-developers/cgat-core>`_ pipelining system. From more information please see the `CGAT-core paper <https://doi.org/10.12688/f1000research.18674.2>`_. Here we illustrate how the pipelines can be run using a toy example which is included in this repository.

Following installation, to find the available pipelines run::

  cellhub -h

Next generate the configuration yml file (for the example pipeline it is empty)::

  cellhub example config -v5

To fully run the example cellhub pipeline run::

  cellhub example make full -v5

However, it may be best to run the individual tasks of the pipeline to get a feel of what each task is doing::

  cellhub example make exampleOriginate -v5

You can also run the pipeline with more advanced combinatorics by running the task::

  cellhub example make advancedRuffus -v5


Getting Started
---------------

To get started please see the :doc:`IFNb example<IFNbExample>`.
