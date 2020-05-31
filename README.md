# cellhub

Tools for the ingest, tracking and retrieval of single-cell data for downstream analysis


## Documentation

When documentation is ready add links here

## Installation

You can manually install cellhub by::

    git clone https://github.com/COMBATOxford/cellhub.git
    cd cellhub
    python setup.py install
    cellhub --help

## Usage

Run the ``cellhub --help`` command to view the help documentation and find available pipelines
to run cellhub.

To generate the configuration file prior to running a pipeline run::

    cellhub example config -v5

To run the main cellhub pipeline run::

    cellhub example make full -v5
