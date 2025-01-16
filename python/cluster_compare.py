import os
import re
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
import logging
import sys
from matplotlib import pyplot as plt



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

#   --source_anndata=%(source_anndata)
#                    --clusterids=%(cluster_ids)s
#                    --ncomp=%(ncomp)s
#                    --outdir=%(outdir)s
#                    --reductiontype=%(reductiontype)s

parser = argparse.ArgumentParser()
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")
parser.add_argument("--clusterids", default="cluster_ids.tsv", type=str,
                    help="path to the cluster ids file")
parser.add_argument("--reduced_dims_name", default="pca", type=str,
                    help="the X_[name] of the reduced dimension")
parser.add_argument("--outdir",default=1, type=str,
                    help="name of the outdir")
parser.add_argument("--ncomps", default=20, type=int,
                    help="Number of reduced components to include in the knn computation")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ######################## Make the dendrogram ############################## #
# ########################################################################### #

# read in the source data
adata = ad.read_h5ad(args.source_anndata)

# , backed='r') no supported see: https://github.com/scverse/scanpy/issues/3199

# add the cluster_ids
clusters = pd.read_csv(args.clusterids, sep="\t")
clusters.index = clusters["barcode_id"]
adata.obs["cluster_id"] = clusters.loc[adata.obs.index, "cluster_id"].astype("category")

# make the dendrogram in the reduced dimensions space
rdims = "X_" + args.reduced_dims_name

sc.tl.dendrogram(adata, 
                 use_rep=rdims,
                 groupby='cluster_id')

# save the plot.
outfile = os.path.join(args.outdir,
                       "cluster.dendrogram.png")

with plt.rc_context():  # Use this to set figure params like size and dpi
    sc.pl.dendrogram(adata, groupby='cluster_id', show=False)
    plt.savefig(outfile, bbox_inches="tight")
                       

L.info("Complete")
