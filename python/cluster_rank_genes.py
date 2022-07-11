# https://scanpy.discourse.group/t/how-to-calculate-the-average-gene-expression-within-each-cluster/126

import os
import re
import argparse
import anndata as ad
import pandas as pd
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
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")
parser.add_argument("--reduced_dims_name", default="pca", type=str,
                    help="the X_[name] of the reduced dimension")
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")
parser.add_argument("--ncomps", default=20, type=int,
                    help="Number of reduced components to include in the knn computation")
parser.add_argument("--method", default="scanpy", type=str,
                    help="scanpy|hnsw|sklearn (scanpy uses pynndescent)")
parser.add_argument("--k", default=20, type=int,
                    help="number of neighbors")
parser.add_argument("--metric", default="euclidean", type=str,
                    help="the distance metric")
parser.add_argument("--threads", default=4, type=int,
                    help="number of threads")
parser.add_argument("--fullspeed", default=False, action="store_true",
                    help="number of threads")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)