# This script filters out a set of cells from a single h5 file
# and saves the result as and anndata

import os
import sys
import anndata as ad
import scanpy as sc
import pandas as pd
import logging
import argparse
import numpy as np
import cellhub.tasks.h5 as h5

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("extract_cells_from_h5.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--barcodes", default=None, type=str,
                    help=('a gzipped file containing the list of barcodes to extract'
                          ' The file should not have a header.')) 
parser.add_argument("--h5", default=None, type=str,
                    help='the h5 file to extract the cells from.')
parser.add_argument("--outfile",default=None, type=str,
                    help="the output file path")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #

L.info("Reading in barcodes")

barcodes = [x for x in pd.read_csv(args.barcodes, header=None)[0].values]

x = h5.read_h5(args.h5)

# use unique identifiers for the index 
# the index is the ensembl gene name, these are not necessarily unique
x.var_names_make_unique()

L.info("Number of barcodes to extract: " + str(len(barcodes)))

common_barcodes = [y for y in barcodes if y in x.obs.index]
     
L.info("Number of barcodes found:" + str(len(common_barcodes)))

if len(barcodes) > len(common_barcodes):

    L.warning("Number of barcodes not found: " + str(len(barcodes) - len(common_barcodes)))

# slice the anndata to the given barcodes    
x = x[common_barcodes]

# save the anndata
x.write_h5ad(args.outfile)

L.info("complete")
