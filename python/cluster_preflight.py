# https://scanpy.discourse.group/t/how-to-calculate-the-average-gene-expression-within-each-cluster/126

import os
import re
import argparse
import anndata as ad
import scanpy as sc
import scipy
import pandas as pd
import numpy as np
import logging
import sys


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--anndata", default="adata.h5ad", type=str,
                    help="path to the anndata object")
parser.add_argument("--reduced_dims_name", default="X_pca", type=str,
                    help="path to the anndata object")
parser.add_argument("--max_reduced_dims", default="10", type=str,
                    help="the maximum number of reduced dims expected")
parser.add_argument("--conserved", action='store_true',
                    help="find conserved markers")   
parser.add_argument("--conserved_factor", default="None", type=str,
                    help="name of the conserved factor")                                           

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ########################### Read in the data ############################## #
# ########################################################################### #

adata = ad.read_h5ad(args.anndata)

L.info("Checking adata.X")
# adata.X contains the scaled data, it should be
# - be dense 
if scipy.sparse.issparse(adata.X):
    raise ValueError("adata.X is sparse. It should contain a dense matrix" 
    " with the scaled values")

L.info('Checking adata.layers["counts"]')
# adata.layers["counts"] contains the counts, it should
# - exist
# - be sparse
if not "counts" in adata.layers.keys():
    raise ValueError("counts missing from source anndata layers")

if not scipy.sparse.issparse(adata.layers["counts"]):
    raise ValueError("counts layer is not sparse")

# sanity check that counts contains integers.
x = adata.layers["counts"].max(axis=1).toarray()
if not np.all([not (i%1) for i in x]):
    raise ValueError("counts layer contains floats")

L.info('Checking adata.layers["log1p"]')
# adata.layers["log1p"] contains the log1p normalised expression values, 
# it should:
# - be sparse
if not "log1p" in adata.layers.keys():
    raise ValueError("log1p missing from source anndata layers")

if not scipy.sparse.issparse(adata.layers["log1p"]):
    raise ValueError("log1p layer is not sparse")

# sanity check that log1p contains floats.
# - note that it is theoretically possible for this check
#   to be incorrect if the max log1p for each cell is
#   an int, or the values have been clipped to ints.
#   ... if encountered in the real world this will be re-written
x = adata.layers["log1p"].max(axis=1).toarray()
if np.all([not (i%1) for i in x]):
    raise ValueError("log1p layer contains integers")

L.info('Checking reduced dimensions exist')
# adata.obsm check the expected key exists
if not "X_" + args.reduced_dims_name in adata.obsm.keys():
    raise ValueError("Reduced dimensions name does not exist in obsm")

L.info('Checking number of reduced dimensions requested')
nrdims = adata.obsm["X_" + args.reduced_dims_name].shape[1]

if int(args.max_reduced_dims) > nrdims:
    message = (str(args.max_reduced_dims) + " dimensions requested, but "
               "only " + str(nrdims) + " dims are present in "
               + "X_" + args.reduced_dims_name)
    raise ValueError(message)

# adata.obs checks

if args.conserved:

    L.info('Checking conserved factor exists in metadata')

    if not args.conserved_factor in adata.obs.columns:
        raise ValueError("The conserved factor was not found in the obs")

L.info("complete")