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
parser.add_argument("--group_means", default="group.means.tsv", type=str,
                    help="mean gene expression per group")
parser.add_argument("--group_sizes", default="None", type=str,
                    help="the sizes of the groups")
parser.add_argument("--group", default="None", type=str,
                    help="The group of interest")
parser.add_argument("--outfile", default="None", type=str,
                    help="The outfile")


args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ########################### Read in the data ############################## #
# ########################################################################### #

group_means = pd.read_csv(args.group_means, sep="\t", index_col=0)
group_means = group_means.sort_index() # unnecessary
group_means.index = group_means.index.astype(str)

group_sizes = pd.read_csv(args.group_sizes, sep="\t", index_col=0)
group_sizes = group_sizes.sort_index() # unnecessary
group_sizes.index = group_sizes.index.astype(str)

# addition of a pseudocount is necessary to avoid Inf values.
pseudo_count = 1

# get means for group of interest
x = group_means.loc[args.group,:] + pseudo_count

# get means for other groups
other_groups = [x for x in group_means.index if x != args.group]
y = group_means.loc[other_groups,:] + pseudo_count

# compute the other mean
group_sizes_other = group_sizes.loc[other_groups, "count"]
group_fractions_other = group_sizes_other / np.sum(group_sizes_other)

z = y.mul(group_fractions_other,0)
other_mean = z.sum(axis=0)

# fold change vs other mean
FC = x/other_mean

# get the minimum fold change
per_group_fold_changes = 1/y.div(x, axis="columns")
min_FC = per_group_fold_changes.min()

results = pd.DataFrame(index=FC.index)
results["gene_id"] = results.index
results["FC"] = FC.loc[results.index]
results["min_FC"] = min_FC.loc[results.index]

results.to_csv(args.outfile, sep="\t", index=False)

# compute the minimum log2 fold change
L.info("complete")