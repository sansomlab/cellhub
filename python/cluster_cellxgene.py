import anndata as ad
import numpy as np
import pandas as pd
import scipy as sc
import sys
import os
import logging
import argparse
from datetime import date
from cellhub.tasks import cellxgene


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("cluster_cellxgene.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")
parser.add_argument("--obs", default="all", type=str,
                    help=('the cols in obs to be included. Either a '
                          'comma seperated list, or "all" '))
parser.add_argument("--umap", default="umap.tsv.gz", type=str,
                    help="a comma seperated list of umaps to add")
parser.add_argument("--umap_facet_x", default="None", type=str,
                    help="a factor to facet the UMAP layout (x axis")
parser.add_argument("--umap_facet_y", default="None", type=str,
                    help="a factor to facet the UMAP layout (x axis")
parser.add_argument("--cluster_paths", default="cluster_ids.tsv.gz", type=str,
                    help="a comma seperated list of cluster_id files to add")
parser.add_argument("--cluster_names", default="umap.tsv.gz", type=str,
                    help="a comma seperated list of the clustering names")
parser.add_argument("--cluster_split", default="None", type=str,
                    help="a factor by which to split the clusters")
parser.add_argument("--adt", default="adt.h5ad", type=str,
                    help='An optional anndata with normalised adt information in the .X')
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################### Functions ################################# #
# ########################################################################### #

# Data format docs
# https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/data-reqs.md

# sanitise the source anndata

# we want to have
# - normalised log1p transformed expression values in .X
# - UMAP layouts
# - some QC vars
# faceting of UMAPs and Clusters by metadata variables of intereset


# Dealing with colors
# - >>> category = "louvain"
# >>> # colors stored in adata.uns must be matplotlib-compatible color information
# >>> adata.uns[f"{category}_colors"]
# array(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22'], dtype='>> # there must be a matching category in adata.obs
# >>> category in adata.obs
# True

L.info("reading the reference anndata object")
source_data = ad.read_h5ad(args.source_anndata, backed='r')

# Add the normalised data
data = ad.AnnData(X=sc.sparse.csc_matrix(source_data.layers["log1p"].copy()))
data.var.index = source_data.var.index.copy()

# Add the metadata from the obs
# populate with columns from cell hub
# and optional factors of interest.
if args.obs == "all":
    data.obs = source_data.obs.copy()
    
    # clean out columns with ":" (duplicates from sql.)
    # remove barcode_id
    data.obs = data.obs[[x for x in data.obs.columns if not ":" in x]]
    
else:
    metadata_cols = [x.strip() for x in args.obs.split(",")]
    metadata_cols = [x for x in metadata_cols if x in data.obs.columns]
    data.obs = source_data.obs[metadata_cols].copy()


L.info("Adding the umap layers")
# Add the umap layers
umap_data = pd.read_csv(args.umap, sep="\t")
umap_data.index = umap_data.barcode_id
data.obsm["X_umap"] = np.array(umap_data.loc[data.obs.index][["UMAP_1","UMAP_2"]])

# Facet the umap
if args.umap_facet_x != "None":
    ux = args.umap_facet_x.strip()
    if ux not in data.obs.columns:
        raise ValueError("umap facet x factor not found in obs")
else:
    ux = None
        
if args.umap_facet_y != "None":
    uy = args.umap_facet_y.strip()
    if uy not in data.obs.columns:
        raise ValueError("umap facet y factor not found in obs")
else:
    uy = None

if ux != None or uy !=None:
    cellxgene.facet_layout(data, 
                 layout="X_umap",
                 name=None,
                 x_factor=ux, 
                 y_factor=uy)


L.info("Adding the clusters")
# Add the clusters 
cluster_paths = [ x.strip() for x in args.cluster_paths.split(",") ]
cluster_names = [ x.strip() for x in args.cluster_names.split(",") ]

if len(cluster_paths) != len(cluster_names):
    raise ValueError("Different number of clusters and cluster names.")

for cluster in zip(cluster_names, cluster_paths):
    L.info("Adding cluster " + cluster[0])
    cluster_data = pd.read_csv(cluster[1], sep="\t")
    cluster_data.index = cluster_data.barcode_id
    
    data.obs[cluster[0]] = cluster_data.loc[data.obs.index]["cluster_id"].astype("category")
    
    L.info("Adding the colors")
    cluster_colors_path = os.path.join(os.path.dirname(cluster[1]),
                                       "cluster_colors.tsv")
    
    if os.path.exists(cluster_colors_path):   
        L.info("reading in the colors")  
        print(cluster_colors_path)
        cluster_colors = pd.read_csv(cluster_colors_path,header=None)[0]
        data.uns[cluster[0] + "_colors"] = np.array(cluster_colors)
    
    if args.cluster_split != "None":
        L.info("splitting the clusters")
        split_clusters = [ "_".join(x) for x in zip(data.obs[cluster[0]].values.astype("str"), 
                                                   data.obs[args.cluster_split].values)]
        split_name = cluster[0] + "_" + args.cluster_split
        data.obs[split_name] = split_clusters
        data.obs[split_name] = data.obs[split_name].astype("category")

L.info("adding the date")
data.uns["ObjectCreateDate"] = date.today().strftime("%Y-%m-%d")

# ########################################################################### #
# ################################### save and exit ######################### #
# ########################################################################### #

L.info("saving the object (with compression")
data.write_h5ad(args.outfile,
                 compression="gzip")

L.info("all done")
