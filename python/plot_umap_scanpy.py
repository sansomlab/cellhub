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
import logging
import yaml

# ########################################################################### #
# ############### Set up the log and figure folder ########################## #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
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

L.info("Running with options ---> %s", opt)

L.info("Writing output to file %s", results_file)

# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #

# read data from h5 file
adata = sc.read_h5ad(results_file)

L.info("Read anndata object")

# ########################################################################### #
# ################################# Run UMAP ################################ #
# ########################################################################### #

if opt["tool"] == 'harmony':
    obsm_use = 'X_harmony'
    key_add = 'harmony'
else:
    obsm_use = 'X_pca'
    key_add = None

dim_red = str(obsm_use.split("_")[1])
L.info("UMAP is run on the following dimension reduction components: " + dim_red)

# 15 neighbors is default, use 30 harmony components or pcas here
sc.pp.neighbors(adata, use_rep=obsm_use , n_pcs = 30, n_neighbors = 15,
                key_added = key_add)
# umap uses the neighbor coordinates 
sc.tl.umap(adata, neighbors_key = key_add)
adata.write(results_file)

umap_coord = pd.DataFrame(adata.obsm['X_umap'])
umap_coord.columns = ["UMAP_1", "UMAP_2"]
umap_coord['barcode'] = adata.obs['barcode'].tolist()
umap_coord.to_csv(os.path.join(opt["outdir"], "umap.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

L.info("Finished running UMAP and writing coordinates")

# ########################################################################### #
# ########################## Make plots ##################################### #
# ########################################################################### #

# split var is part of plot_vars so at least one variable to plot
L.info("Plot variables on UMAP")

if ',' in opt["plot_vars"]:
    for v in opt["plot_vars"].split(','):
        L.info("Making plot for variable: " + str(v))
        file_name = "_" + str(v)
        sc.pl.umap(adata, color=str(v), save = file_name + ".png", show=False)
        sc.pl.umap(adata, color=str(v), save = file_name + ".pdf", show=False)
else:
    L.info("Making plot for variable: " + str(opt["plot_vars"]))
    file_name = "_" + str(opt["plot_vars"])
    sc.pl.umap(adata, color=str(opt["plot_vars"]), save = file_name + ".png", show=False)
    sc.pl.umap(adata, color=str(opt["plot_vars"]), save = file_name + ".pdf", show=False)

L.info("Done UMAP plotting")

L.info("Completed")
