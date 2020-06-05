"""===========================
Pipeline fetch cells
===========================

:Author: Sansom lab
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
========

This pipeline fetches a given set of cells from market matrices or
loom files into a single market matrix file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_fetch_cells.py config


Input files
-----------

There a two inputs, both specified in the pipeline.yml file.

(1) A map of cell barcodes to matrix identifiers

This file should be a gzipped tsv file with two columns. Column
1 should have the title "barcode" and contain the cell barcodes.
Column 2 should have the title "matrix_id" and contain a matrix
identifier

barcode     matrix_id
ACCCATCG    channel_1
ATTCATCG    channel_1
AGGCATCG    channel_2
TCCCATCG    channel_2
gCCCATCG    channel_3


(2) A table containing the matrix identifiers, locations and types.

This file should be a tsv file with three columns:

matrix_id    matrix_dir                    matrix_type
channel_1    /full/path/channel_1_matrix   mm
channel_2    /full/path/channel_2_matrix   mm
channel_3    /full/path/channel_3_matrix   mm

The matrix type should be specified as either "mm" for market
matrix format or "loom" for the loom file format.


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
import hashlib
from pathlib import Path
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P


# -------------------------- < parse parameters > --------------------------- #

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


# set the location of the tenx code directory
if "cellhub_dir" not in PARAMS.keys():
    PARAMS["cellhub_dir"] = Path(__file__).parents[1]
else:
    raise ValueError("Could not set the location of the "
                     "cellhub code directory")


# ----------------------- < pipeline configuration > ------------------------ #

# handle pipeline configuration
if len(sys.argv) > 1:
        if(sys.argv[1] == "config") and __name__ == "__main__":
                    sys.exit(P.main(sys.argv))


# ------------------------------ < functions > ------------------------------ #

# from:
# https://stackoverflow.com/questions/3431825/
# generating-an-md5-checksum-of-a-file

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# ########################################################################### #
# ############## read in the cell and matrix information #################### #
# ########################################################################### #

# get the map of cells to matrix identifies
cell_map = PARAMS["cell_to_matrix_map_file"]
cells = pd.read_csv(cell_map, sep="\t")

# get the table containing the matrix information
matrix_table = PARAMS["matrix_table_file"]
matrices = pd.read_csv(matrix_table, sep="\t")
matrices.index = matrices["matrix_id"]


# ########################################################################### #
# ############################# pipeline tasks ############################## #
# ########################################################################### #

def matrix_subset_jobs():
    '''Generate a list of subsets from the cell manifest'''

    matrix_ids = [x for x in cells.matrix_id.unique()]

    for matrix_id in matrix_ids:

        matrix_subset_dir = os.path.join("matrix.subsets.dir",
                                         matrix_id)

        matrix_subset_barcodes = os.path.join(matrix_subset_dir,
                                              "cells_to_extract.txt.gz")

        yield [None, matrix_subset_barcodes]


@follows(mkdir("matrix.subsets.dir"))
@files(matrix_subset_jobs)
def setupSubsetJobs(infile, outfile):
    '''Setup the folders for the subsetting jobs'''

    matrix_subset_dir = os.path.dirname(outfile)
    matrix_id = os.path.basename(Path(outfile).parent)

    os.mkdir(matrix_subset_dir)

    cell_subset = cells["barcode"][cells["matrix_id"] == matrix_id]

    cell_subset.to_csv(outfile,
                       header=False,
                       index=False,
                       compression="gzip")


@transform(setupSubsetJobs,
           regex(r"matrix.subsets.dir/(.*)/cells_to_extract.txt.gz"),
           r"matrix.subsets.dir/\1/matrix.mtx.gz")
def cellSubsets(infile, outfile):
    '''Extract a given subset of cells from a matrix'''

    to_cluster = True

    # matrix_subset_dir = os.path.dirname(outfile)
    matrix_id = os.path.basename(Path(outfile).parent)

    outdir = os.path.dirname(infile)

    matrix_dir = matrices.loc[matrix_id]["matrix_dir"]

    matrix_type = matrices.loc[matrix_id]["matrix_type"]

    statement = '''Rscript %(cellhub_dir)s/R/extract_cells.R
                           --cells=%(infile)s
                           --matrixdir=%(matrix_dir)s
                           --matrixtype=%(matrix_type)s
                           --outdir=%(outdir)s
                '''

    P.run(statement)


@follows(mkdir("output.dir"))
@merge(cellSubsets,
       "output.dir/matrix.mtx.gz")
def mergeSubsets(infiles, outfile):
    '''merge the market matrix files into a single matrix'''

    ncells = 0

    # get the dimensions of all the market matrix files
    mtx_specs = {}

    mtx_encoding = "us-ascii"

    for subset in infiles:

        matrix_id = os.path.basename(Path(subset).parent)

        with gzip.open(subset, "r") as mtx:
            for i, line in enumerate(mtx, 1):
                if i == 1:
                    mtx_header = line  # .decode("us-ascii").strip()
                if i == 2:
                    line_str = line.decode(mtx_encoding)
                    nrow, ncol, nnonzero = line_str.strip().split(" ")
                    mtx_specs[matrix_id] = {"nrow": int(nrow),
                                            "ncol": int(ncol),
                                            "nnonzero": int(nnonzero),
                                            "mtx_file": subset}
                    break

    matrices_to_merge = tuple(mtx_specs.keys())

    mtx_outfile = outfile[:-len(".gz")]
    barcodes_outfile = os.path.join(os.path.dirname(outfile),
                                    "barcodes.tsv")

    features_outfile = os.path.join(os.path.dirname(outfile),
                                    "features.tsv.gz")

    if os.path.exists(mtx_outfile) or \
       os.path.exists(mtx_outfile + ".gz"):
        raise ValueError("mtx outfile already exists")

    if os.path.exists(barcodes_outfile) or \
       os.path.exists(barcodes_outfile + ".gz"):
        raise ValueError("barcodes outfile already exists")

    if os.path.exists(features_outfile):
        raise ValueError("features outfile already exists")

    # construct the header of the market matrix file.
    nrows = []
    total_ncol = 0
    total_nnonzero = 0

    for matrix_id in matrices_to_merge:
        nrows.append(mtx_specs[matrix_id]["nrow"])
        total_ncol += mtx_specs[matrix_id]["ncol"]
        total_nnonzero += mtx_specs[matrix_id]["nnonzero"]

    if len(set(nrows)) > 1:
        raise ValueError("the input matrices have different numbers of rows!")

    out_spec = " ".join([str(nrows[0]),
                         str(total_ncol),
                         str(total_nnonzero)]) + "\n"

    # write the header of the market matrix file
    with open(mtx_outfile, "wb") as out:
        out.write(mtx_header)
        out.write(out_spec.encode(mtx_encoding))

    column_offset = 0

    feature_file_checksums = []

    for matrix_id in matrices_to_merge:

        mtx_file = mtx_specs[matrix_id]["mtx_file"]

        barcodes_file = os.path.join(os.path.dirname(mtx_file),
                                     "barcodes.tsv.gz")

        # append the matrix values, offsetting the column index
        statement = '''zcat %(mtx_file)s
                       | awk 'NR>2{print $1,$2 + %(column_offset)s,$3}'
                       >> %(mtx_outfile)s
                    '''

        P.run(statement)

        # append the barcodes, adding the matrix identifier
        statement = '''zcat %(barcodes_file)s
                       | awk '{print $1"-%(matrix_id)s"}'
                       >> %(barcodes_outfile)s
                    '''

        P.run(statement)

        # increase the column offset by the number of appended columns
        column_offset += mtx_specs[matrix_id]["ncol"]

        features_file = os.path.join(os.path.dirname(mtx_file),
                                     "features.tsv.gz")

        feature_file_checksums.append(md5(features_file))

    # check that all of the subsets have identical features.
    if len(set(feature_file_checksums)) != 1:
        raise ValueError("The matrices have different features")

    # compress the outfiles
    statement = '''gzip %(mtx_outfile)s;
                   gzip %(barcodes_outfile)s
                '''

    P.run(statement)

    # copy over the features to the output directory
    copyfile(features_file, features_outfile)


# ########################################################################### #
# ##################### full target: to run all tasks ####################### #
# ########################################################################### #

@follows(mergeSubsets)
def full():
    pass


# ------------------- < ***** end of pipeline **** > ------------------------ #

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
