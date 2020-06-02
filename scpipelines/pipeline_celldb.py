'''
===============
Pipeline celldb
===============


Overview
========

This pipeline uploads the outputs from the upstream single-cell preprocessing
steps into a SQLite database.

Usage
=====

See :ref:`` and :ref:`` for general information on how to use cgat
pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

Default configuration files can be generated by executing:
   cellhub celldb config

Input files
-----------

The pipeline requires the output of ??. Each task of the upstream
pipeline generates a tsv configured file.

Dependencies
------------

Pipeline output
===============

The pipeline returns an SQLite populated database of metadata and
quality features that aid the selection of 'good' cells from 'bad' cells.

Currently the following tables are generated:
* metadata


Code
====

'''

from ruffus import *

import sys
import os
import re
import pandas

import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools


# Load options from the config file

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Insert database functions here

# Example of loading all tsv to database
# Need to limit jobs to 1 because cant have concurrent connections with SQLite
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(,
           suffix(".tsv"),
           ".load")
def loadExample(infile, outfile):
    '''load comparison data into database.'''

    # csvdb otpions
    options = ""
    P.load(infile, outfile, options)


# -------------------------------


def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
