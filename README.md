# cellhub

Tools for the ingest, tracking and retrieval of single-cell data for downstream analysis


## Documentation

When documentation is ready add links here

## Installation

You can manually install cellhub by::

    git clone https://github.com/COMBATOxford/cellhub.git
    cd cellhub
    python setup.py develop
    cellhub --help

## Usage

Run the ``cellhub --help`` command to view the help documentation and find available pipelines
to run cellhub.

I have included a example pipeline with a set of ruffus decorators that
demonstrates the functionality of cgatcore pipelines.

Following installation, to find the available pipelines run

    cellhub -h

Next generate the condifuration yml file (for the example pipleine it is empty)

    cellhub example config -v5

To fully run the example cellhub pipeline run

    cellhub example make full -v5

However, it may be best to run the individual tasks of the pipeline to get
a feel of what each task is doing

    cellhub example make exampleOriginate -v5

You can also run the pipeline with more advanced combinatorics
by running the task

    cellhub example make advancedRuffus -v5
