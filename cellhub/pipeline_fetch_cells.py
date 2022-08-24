"""
====================
Pipeline fetch cells
====================

Overview
========

This pipeline fetches a given set of cells from market matrices or
loom files into a single market matrix file.

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` on general
information how to use CGAT pipelines.

Configuration
-------------

It is recommended to fetch the cells into a new directory. Fetching
of multiple datasets per-directory is (deliberately) not supported.

The pipeline requires a configured :file:`pipeline_fetch_cells.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_fetch_cells.py config


Inputs
------

The pipeline will fetch cells from a cellhub instance according to
the parameters specified in the local pipeline_fetch_cell.yml file.


The location of the cellhub instances must be specificed in the yml: ::

   cellhub:
       location: /path/to/cellhub/instance

The specifications of the cells to retrieve must be provided as an SQL
statement (query) that will be executed against the "final" table of the cellhub
database: ::

    cellhub:
        sql_query: >-
            select * from final
            where pct_mitochondrial < 10
            and ngenes > 200;


The cells will then be automatically retrieved from the API.

Dependencies
------------

This pipeline requires:


Pipeline output
===============

The pipeline outputs a folder containing a single market matrix
that contains the requested cells.

"""

import os
import sys
import gzip
from shutil import copyfile

from pathlib import Path
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools

import cellhub.tasks.parameters as chparam
import cellhub.tasks.fetch_cells as fetch_cells

# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = chparam.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(chparam.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]

# ------------------------------ Pipeline Tasks ------------------------------ #


@follows(mkdir("fetch.cells.dir"))
@files(os.path.join(PARAMS["cellhub_location"],"celldb.dir/csvdb"),
       "fetch.cells.dir/cell.table.sentinel")
def fetchCells(infile, outfile):
    '''Fetch the table of the user's desired cells from the database
       effectively, cell-metadata tsv table.
    '''

    cell_table = outfile.replace(".sentinel", ".tsv.gz")

    if(PARAMS["sample"]=="all"):
        sample = ""
    else:
        take = int(PARAMS["sample"]) + 1
        sample = '''| body shuf | awk 'NR <= %(take)s' ''' % locals()

    statement ='''body() {
                   IFS= read -r header;
                   printf '%%s\\n' "$header";
                  "$@";
                  };

              sqlite3 -header %(infile)s
                      -separator $'\\t'
                      '%(cellhub_sql_query)s'
                  %(sample)s
                  | gzip -c
                  > %(cell_table)s'''

    P.run(statement)
    IOTools.touch_file(outfile)


# ########################################################################### #
# ######################### fetch GEX data ################################## #
# ########################################################################### #


@follows(mkdir("anndata.dir"))
@files(fetchCells,
       "anndata.dir/gex.sentinel")
def GEX(infile, outfile):
    '''
    Extract the target cells into a single anndata.
    Note that this currently contains all the modalities

    TODO: support down-sampling
    '''

    cell_table = infile.replace(".sentinel", ".tsv.gz")
    api_path = os.path.join(PARAMS["cellhub_location"],"api")
    outdir = os.path.dirname(outfile)   
    log_file = outfile.replace(".sentinel", ".log")

    statement = '''python %(cellhub_code_dir)s/python/fetch_cells_from_h5.py 
                   --cells=%(cell_table)s
                   --feature_type=GEX
                   --api=%(api_path)s
                   --outname=gex.h5ad
                   --outdir=%(outdir)s
                   &> %(log_file)s
                '''

    P.run(statement)
    IOTools.touch_file(outfile)


# # ########################################################################### #
# # ######################### fetch ADT data ################################## #
# # ########################################################################### #

@follows(mkdir("anndata.dir"))
@active_if(PARAMS["ADT_fetch"])
@files(fetchCells,
       "anndata.dir/adt.sentinel")
def ADT(infile, outfile):
    '''
    Extract the target cells into a single anndata.
    Note that this currently contains all the modalities

    TODO: support down-sampling
    '''

    cell_table = infile.replace(".sentinel", ".tsv.gz")
    api_path = os.path.join(PARAMS["cellhub_location"],"api")
    outdir = os.path.dirname(outfile)   
    log_file = outfile.replace(".sentinel", ".log")

    statement = '''python %(cellhub_code_dir)s/python/fetch_cells_from_h5.py 
                   --cells=%(cell_table)s
                   --feature_type=ADT
                   --api=%(api_path)s
                   --outdir=%(outdir)s
                   --outname=adt.h5ad
                   &> %(log_file)s
                '''

    P.run(statement)
    IOTools.touch_file(outfile)


# ########################################################################### #
# ######################### fetch VDJ_B data ################################ #
# ########################################################################### #

# Fetch BCR data from the api here.


# ########################################################################### #
# ######################### fetch VDJ_T data ################################ #
# ########################################################################### #

# Fetch TCR data from the api here.


# ########################################################################### #
# ##################### full target: to run all tasks ####################### #
# ########################################################################### #


@follows(GEX, ADT)
def full():
    pass


# ------------------- < ***** end of pipeline **** > ------------------------ #


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
