import os
import re
import argparse
import anndata as ad
#import scanpy as sc
#import scipy
#import pandas as pd
#import numpy as np
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
parser.add_argument("--loom", default="X_pca", type=str,
                    help="path to save the loom object")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ########################### Read in the data ############################## #
# ########################################################################### #

adata = ad.read_h5ad(args.anndata)

if(len(set(adata.var.index.values))) < len(adata.var.index.values):
    L.info("Making var names unique")
    adata.var_names_make_unique()

# the gene_name is taken as the var.index
adata.var["gene_name"] = adata.var.index.copy()

adata.write_loom(args.loom)

L.info("complete")