import os
import sys
import scanpy
import pandas as pd
import logging
import argparse
import numpy as np


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("extract_pseudobulks.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--anndata", default=None, type=str,
                    help="a h5ad anndata object")
parser.add_argument("--layer", default=".X", type=str,
                    help='location of the count data, ".X" or layer name')
parser.add_argument("--clusters",default=None, type=str,
                    help="a tsv file with columns barcode and cluster_id")
parser.add_argument("--samplecol",default="sample", type=str,
                    help="Name of a column in the .obs containing the sample names")
parser.add_argument("--manifest",default=None, type=str,
                    help='An optional tsv mapping of "cluster_id" to "group" names')
parser.add_argument("--name", default="pseudocounts", type=str,
                    help="a name for the result")
parser.add_argument("--outdir", default=".", type=str,
                    help="the directory for the output")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #


## Read in the adata

x = scanpy.read_h5ad(args.anndata)

## read in the clusters

clust = pd.read_csv(args.clusters, sep="\t")

clust.index = clust["barcodes"]

## read in the manifest

manifest = pd.read_csv(args.manifest, sep="\t")

samples = [ x for x in x.obs[args.samplecol].unique()]

groups = [x for x in manifest["group"].unique()]

ncol = len(samples) * len(groups)
nrow = len(x.var.index)

pseudocounts = np.zeros((nrow, ncol))

col_idx = 0
col_names = []

for group in groups:

    L.info("Working on group: " + group)

    clusters = manifest["cluster_id"][manifest["group"]==group].values

    # L.info("Fetching barcodes for clusters: " + ",".join([str(x) for x in clusters]))

    barcodes_to_fetch = clust.index[clust["cluster_id"].isin( clusters)]

    # L.info("Found " + str(len(barcodes_to_fetch)) + " barcodes")

    for sample in samples:

           # L.info("Working on sample: " + sample)
           sample_barcodes = x.obs.index[x.obs.index.isin(barcodes_to_fetch) &
                                         (x.obs[args.samplecol]==sample)]

           # L.info("Got " + str(len(sample_barcodes)) + " sample barcodes for this group")

           if args.layer == ".X":
               # L.info("Summing on matrix in .X")
               pseudocounts[:, col_idx] = np.sum(x[sample_barcodes].X,0)

           else:
               # L.info("Summing in matrix in layer: " + args.layer)
               pseudocounts[:, col_idx] = np.sum(x[sample_barcodes].layers[args.layer],0)

           col_names.append(str(group) + "_" + str(sample))
           col_idx += 1

L.info("Summation of counts complete")
pseudocount_frame = pd.DataFrame(pseudocounts,
                                 columns = col_names,
                                 index = x.var.index)

pseudocount_frame = pseudocount_frame[col_names].astype(int)

pseudocount_frame["gene"] = pseudocount_frame.index

pseudocount_frame = pseudocount_frame[["gene"] + col_names]

pseudocount_frame.to_csv(os.path.join(args.outdir, args.name)
                         + ".tsv.gz", sep="\t",
                         index=False, header=True)

L.info("complete")
