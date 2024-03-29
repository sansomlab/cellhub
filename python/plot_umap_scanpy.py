#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import os
import scanpy as sc
import harmonypy as hm
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import logging
import yaml
warnings.filterwarnings('ignore')


def removeMetadata(adata, metadata_infile, id_col):
    metadata_cols = pd.read_csv(metadata_infile, sep = "\t",
                                          compression="gzip")
    cols_remove = metadata_cols.columns.copy()
    cols_remove = cols_remove.drop(labels=[id_col])
    adata.obs.drop(cols_remove, inplace=True, axis = 1)
    return adata

# ########################################################################### #
# ############### Set up the log and figure folder ########################## #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.WARNING)
L = logging.getLogger("run_umap")

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# comment this because of numba issues
#sc.logging.print_versions()

# ########################################################################### #
# ############################ Script arguments ############################# #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--task-yml", default="", type=str,
                    help="yml for this task")
args = parser.parse_args()


# Read YAML file
with open(args.task_yml, 'r') as stream:
    opt = yaml.safe_load(stream)

# figures folder
sc.settings.figdir = os.path.join(opt["outdir"], "figures.dir")
sc.settings.set_figure_params(dpi=300, dpi_save=300)

# write folder
results_file = os.path.join(opt["outdir"], "normalized_integrated_anndata.h5ad")

L.warning("Running with options ---> %s", opt)

L.warning("Writing output to file %s", results_file)

# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #

# read data from h5 file
adata = sc.read_h5ad(results_file)

L.warning("Read anndata object")

# ########################################################################### #
# ################################# Run UMAP ################################ #
# ########################################################################### #

if opt["tool"] == 'harmony':
    obsm_use = 'X_harmony'
    key_add = 'harmony'
elif opt["tool"] == "bbknn":
    obsm_use = 'X_pca'
    key_add = None
elif opt["tool"] == "scanorama":
    obsm_use = 'X_scanorama_embedding'
    key_add = 'scanorama'
else:
    obsm_use = 'X_pca'
    key_add = None

dim_red = str(obsm_use.split("_")[1])
L.warning("UMAP is run on the following dimension reduction components: " + dim_red)

# 15 neighbors is default, use 30 harmony components or pcas here
if opt["tool"] == "bbknn":
    L.warning("No need to determine neighbors again")
else:
    L.warning("Find neighbors")
    sc.pp.neighbors(adata, use_rep=obsm_use, n_neighbors = 15,
                    key_added = key_add)
# umap uses the neighbor coordinates
sc.tl.umap(adata, neighbors_key = key_add, min_dist=opt["umap_min_dist"] )
adata.write(results_file)
L.warning("Finished running UMAP")


umap_coord = pd.DataFrame(adata.obsm['X_umap'])
umap_coord.columns = ["UMAP_1", "UMAP_2"]
umap_coord['barcode'] = adata.obs.index.tolist()
umap_coord.to_csv(os.path.join(opt["outdir"], "umap.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

L.warning("Finished and writing UMAP coordinates")

# write out metadata for plotting in R
metadata_out = adata.obs.copy()
# if barcode exists, this is likely incorrect -> use index instead
if 'barcode' in metadata_out.columns:
    del metadata_out['barcode']
metadata_out.index.name = 'barcode'
metadata_out.reset_index(inplace=True)
metadata_out.head()
metadata_out.to_csv(os.path.join(opt["outdir"], "metadata.tsv.gz"),
                    sep="\t", index=False, compression="gzip")


# ########################################################################### #
# ########################## Make plots ##################################### #
# ########################################################################### #

# split var is part of plot_vars so at least one variable to plot
L.warning("Plot variables on UMAP")

if ',' in opt["plot_vars"]:
    for v in opt["plot_vars"].split(','):
        L.warning("Making plot for variable: " + str(v))
        file_name = "_" + str(v)
        if v not in adata.obs.columns:
            L.warning("Variable is not in the metadata.")
        else:
            sc.pl.umap(adata, color=str(v), save = file_name + ".png", show=False)
            #sc.pl.umap(adata, color=str(v), save = file_name + ".pdf", show=False)
else:
    L.warning("Making plot for variable: " + str(opt["plot_vars"]))
    file_name = "_" + str(opt["plot_vars"])
    if opt["plot_vars"] not in adata.obs.columns:
        L.warning("Variable is not in the metadata.")
    else:
        sc.pl.umap(adata, color=str(opt["plot_vars"]), save = file_name + ".png", show=False)
        #sc.pl.umap(adata, color=str(opt["plot_vars"]), save = file_name + ".pdf", show=False)

L.warning("Done UMAP plotting")

# remove metadata from anndata before writing
if 'metadata_file' in opt.keys():
    adata = removeMetadata(adata = adata,
                           metadata_infile = opt["metadata_file"],
                           id_col = opt["metadata_id"])

adata.write(results_file)
L.warning("Removed metadata columns")

L.warning("Completed")
