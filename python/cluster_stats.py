# https://scanpy.discourse.group/t/how-to-calculate-the-average-gene-expression-within-each-cluster/126

import os
import re
import argparse
import anndata as ad
import scanpy as sc
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
parser.add_argument("--subset_factor", default="None", type=str,
                    help="A column name in obs for subsetting the cells")
parser.add_argument("--subset_level", default="None", type=str,
                    help="The level to use to extract the cell subset")
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")
parser.add_argument("--clusterids", default="cluster_ids.tsv.gz", type=str,
                    help="cluster identifiers")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ########################### Read in the data ############################## #
# ########################################################################### #

source_adata = ad.read_h5ad(args.anndata, backed='r')

L.info("Performing expm1 computation")
# we want to do this on the untransformed total count normalised expression values
adata = ad.AnnData(X=source_adata.layers["log1p"].expm1().copy(),
                obs=source_adata.obs.copy(), var=source_adata.var.copy())

source_adata.file.close()

cids = pd.read_csv(args.clusterids, sep="\t")
cids.index = cids["barcode_id"]

adata.obs["cluster_id"] = cids.loc[adata.obs.index, "cluster_id"].astype("category")

if not args.subset_factor == "None":
    L.info('subsetting by factor "' +
    args.subset_factor + '" to level "' + args.subset_level + '"')
    
    L.info('data shape before subsetting')
    print(adata.X.shape)
    adata = adata[adata.obs[args.subset_factor]==args.subset_level].copy()
    
    L.info('data shape after subsetting')
    print(adata.X.shape)

# ########################################################################### #
# ######################### Compute statistics ############################## #
# ########################################################################### #

L.info("computing cluster means")
res = pd.DataFrame(columns=adata.var_names, index=adata.obs['cluster_id'].cat.categories)                                                                                                 
for clust in adata.obs.cluster_id.cat.categories: 
    res.loc[clust] = adata[adata.obs['cluster_id'].isin([clust]),:].X.mean(0)

L.info("saving cluster means")
res.to_csv(args.outfile, sep="\t")

L.info("saving cluster sizes")
cell_number_file = args.outfile.replace("stats.tsv.gz","sizes.tsv.gz")
cell_numbers = pd.DataFrame(adata.obs.cluster_id.value_counts())
cell_numbers.columns = ["count"]
cell_numbers.to_csv(cell_number_file, sep="\t")

L.info("complete")