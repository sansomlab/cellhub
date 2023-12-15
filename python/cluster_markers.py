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
parser.add_argument("--test_factor", default="None", type=str, 
                    help="A column name in obs for within cluster a vs b tests")
parser.add_argument("--a", default="None", type=str,
                    help="The level name of group a") 
parser.add_argument("--b", default="None", type=str,
                    help="The level name of group a")
parser.add_argument("--group_means", default="group.means.tsv", type=str,
                    help="mean gene expression per group")
parser.add_argument("--group_sizes", default="None", type=str,
                    help="the sizes of the groups")    
parser.add_argument("--pseudocount", default=0.1, type=float,
                    help="pseudocount used for fold change estimation")                          
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")
parser.add_argument("--clusterids", default="cluster_ids.tsv.gz", type=str,
                    help="cluster identifiers")
parser.add_argument("--cluster", default="0", type=str,
                    help="the cluster to find markers for")
parser.add_argument("--method", default="t-test", type=str,
                    help="scanpy rank_genes_groups method")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# ########################################################################### #
# ########################### Read in the data ############################## #
# ########################################################################### #

L.info("Reading in the source anndata:")
source_adata = ad.read_h5ad(args.anndata, backed='r')

adata = ad.AnnData(var = source_adata.var, obs = source_adata.obs)
adata.layers["counts"] = source_adata.layers["counts"]
adata.layers["log1p"] = source_adata.layers["log1p"]

cids = pd.read_csv(args.clusterids, sep="\t", dtype=str)
cids.index = cids["barcode_id"]
cids["cluster_id"] = cids["cluster_id"].astype("category")
adata.obs["cluster_id"] = cids.loc[adata.obs.index, "cluster_id"]

if not args.subset_factor == "None":
    L.info('subsetting by factor "' +
    args.subset_factor + '" to level "' + args.subset_level + '"')
    L.info('data shape before:')
    print(adata.layers["counts"].shape)
    adata = adata[adata.obs[args.subset_factor]==args.subset_level].copy()
    L.info('data shape after:')
    print(adata.layers["counts"].shape)


# ########################################################################### #
# ######################### Compute statistics ############################## #
# ########################################################################### #

# scanpy.tl.rank_genes_groups(adata, groupby, use_raw=None, groups='all', 
# reference='rest', n_genes=None, rankby_abs=False, pts=False, 
# key_added=None, copy=False, method=None, corr_method='benjamini-hochberg', 
# tie_correct=False, layer=None, **kwds)

# deal with case where base key is missing
# in the adata.uns['log1p'] dictionary.
if 'log1p' in adata.uns_keys():
    if 'base' not in adata.uns['log1p'].keys():
        adata.uns['log1p']['base'] = None

# cluster vs other test.
if args.test_factor == "None":
    L.info("discovering markers for cluster " + str(args.cluster))
    
    L.info("Checking if cluster is present in data:") #this is important for conserved markers
    #check if more than 5 cells
    if(adata[adata.obs['cluster_id'].isin([args.cluster]),].shape[0] > 5):
    
        sc.tl.rank_genes_groups(adata,
                                use_raw=False,
                                layer="log1p",
                                groupby="cluster_id",
                                groups = [args.cluster],
                                method = args.method,
                                pts = True)
    else:
        L.info("Cluster has less than 5 cells. No markers discovered")


else:
    L.info("performing a vs b test on factor " + args.test_factor)
    # subset the cells to the cluster of interest
    adata = adata[adata.obs["cluster_id"]==str(args.cluster)].copy()
    # subset to only the two groups of interest
    # Strangely, this is necessary for computation of pct in the other group.
    adata = adata[adata.obs[args.test.factor].isin([args.a, args.b])].copy()

    sc.tl.rank_genes_groups(adata,
                            use_raw=False,
                            layer="log1p",
                            groupby= args.test_factor,
                            groups = [args.a],
                            method = args.method,
                            pts = True)


if(adata[adata.obs['cluster_id'].isin([args.cluster]),].shape[0] > 5):
    # save the result.
    L.info("formatting results matrix")

    xx = adata.uns["rank_genes_groups"]

    # Note: the pts columns have their own index and a different row order!
    #       here we need to work around this sub-optimal situation.

    res = pd.DataFrame({"gene_name":[x[0] for x in xx["names"]]})
    res.index = res["gene_name"]

    if "gene_ids" in adata.var.columns:
        # Add the gene_ids to the result
        res["gene_id"] = adata.var.loc[res["gene_name"],"gene_ids"]

    res["pval"] = [x[0] for x in xx["pvals"]]
    res["score"] = [x[0] for x in xx["scores"]]
    res["logfoldchange"] = [x[0] for x in xx["logfoldchanges"]]
    res["pct"] = xx["pts"].loc[res.index, args.cluster].values
    res["pct_other"] = xx["pts_rest"].loc[res.index, args.cluster].values


    # ########################################################################### #
    # ################### compute the fold changes ############################## #
    # ########################################################################### #

    # scanpy does this:
    #   foldchanges = (self.expm1_func(mean_group) + 1e-9) / (
    #                 self.expm1_func(mean_rest) + 1e-9
    #             )  # add small value to remove 0's

    # https://github.com/scverse/scanpy/blob/1c9e404b7462ee52a0664517c0f569ff16791ca4/scanpy/tools/_rank_genes_groups.py#L417-L419
    #
    # Seurat does this to calculate the mean expression of each group ("data" assay):
    #  log(x = rowMeans(x = expm1(x = x)) + pseudocount.use   # default pseudocount.use = 1
    #  
    # Some things are obvious:
    # (1) Scanpy authors are prioritising speed over accuracy by doing expm(mean) 
    #     instead of mean(expm). This doesn't seem necessary.
    # (2) The Scanpy pseudocount is extremely small. A cluster of 100 cells with a 
    #     gene expressed in 1 cell will get a raw fold change of 1e+07 (10 million) against
    #     clusters with no expression!
    #
    # What should the pseudocount be set to?
    # - Considerations
    #    - 10x polyA capture is perhaps 5-10% efficient
    #    - So if 5% of cells in cluster have a marker that is not expressed
    #      elsewhere, we would want to report this marker (assuming BH p.adj > 0.5)
    #    - (5/100 + 1e-9) / 1e-9 = 5e7 # scapy fold change - doesn't make biological sense.
    #    - (5/100 + 1) / 1 = 1.05 # seurat fold change - real signal might be missed?
    #    - (5/100 + 0.1) / 0.1 = 1.5 # seems reasonable
    #    
    # Based on this thought experiment Seurat seems perhaps a little conservative (although
    # we do not consider noise) so we recommend to use 
    # a pseudocount of 0.1 and require a minimum raw fold change of >=1.5. 
    #
    # Here we estimate the true fold changes on the mean(expm1) values.
    # - these are are expected to be similar but not exactely the same
    #   as the scanpy "approximations"!
    # (quick arbitrary check of sig. markers gives pearson r 0.998, spearmans rho 0.994)

    L.info("Reading the group means")
    group_means = pd.read_csv(args.group_means, sep="\t", index_col=0)
    group_means = group_means.sort_index() # unnecessary
    group_means.index = group_means.index.astype(str)

    L.info("Reading the group sizes")
    group_sizes = pd.read_csv(args.group_sizes, sep="\t", index_col=0)
    group_sizes = group_sizes.sort_index() # unnecessary
    group_sizes.index = group_sizes.index.astype(str)

    L.info("Adding the pseudocount")
    # addition of a pseudocount is necessary to avoid Inf values.
    # as above it is important to get this right if we wish to 
    # preserve some biological intuition. 
    pseudo_count = args.pseudocount  # scanpy 1e-9; seurat = 1.

    # get means for group of interest
    x = group_means.loc[args.cluster,:] + pseudo_count

    L.info("Computing the mean of the other group")
    # get means for other groups
    other_groups = [x for x in group_means.index if x != args.cluster]
    y = group_means.loc[other_groups,:] 

    # compute the other mean
    group_sizes_other = group_sizes.loc[other_groups, "count"]
    group_fractions_other = group_sizes_other / np.sum(group_sizes_other)

    z = y.mul(group_fractions_other,0)
    other_mean = z.sum(axis=0)

    other_mean = other_mean + pseudo_count

    L.info("Computing the fold changes")
    print(x.shape)
    print(other_mean.shape)

    # fold change vs other mean
    FC = x/other_mean

    print(FC)

    # get the minimum fold change
    per_group_fold_changes = 1/y.div(x, axis="columns")
    min_FC = per_group_fold_changes.min()

    #results = pd.DataFrame(index=FC.index)
    #results["gene_id"] = results.index

    # (the mean column cannot be called "mean" due to issues with dplyr)
    res["mean_exprs"] = x - pseudo_count
    res["mean_exprs_other"] = other_mean - pseudo_count
    res["FC"] = FC.loc[res.index]
    res["min_FC"] = min_FC.loc[res.index]

    L.info("saving the results")
    res.to_csv(args.outfile, sep="\t",index=False)

    L.info("complete")