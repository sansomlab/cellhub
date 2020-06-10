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
L = logging.getLogger("run_harmony")

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
sc.settings.figdir = opt["outdir"]
sc.settings.set_figure_params(dpi=300, dpi_save=300)

# write folder                                                                                    
results_file = os.path.join(opt["outdir"], "normalized_integrated_anndata.h5ad")

L.info("Running with options ---> %s", opt)

L.info("Writing output to file %s", results_file)

# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #


adata = sc.read_10x_mtx(opt["matrixdir"],
                        var_names='gene_symbols', cache=True)
# could use gene_symbols here but then also required:
adata.var_names_make_unique()

# ## add metadata for cells
metadata = pd.read_csv(os.path.join(opt["matrixdir"], "metadata.tsv.gz"), sep="\t")
adata.obs = metadata

L.info("Made anndata object")

# ########################################################################### #
# ################ Run normalization and scaling ############################ #
# ########################################################################### #

## Basic normalization & log-transformation

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

L.info("Finished normalization and log-transform")

## Identify highly variable genes

# If n_top_genes is used then the flavour = 'cellranger'
if 'hv_genes' in opt.keys():
    L.info("Use hv genes from list")
    gene_list = pd.read_csv(opt["hv_genes"])
    #adata.var.highly_variable = 
else:
    L.info("Determine hv genes using scanpy")
    sc.pp.highly_variable_genes(adata, n_top_genes=opt["ngenes"])
    #sc.pl.highly_variable_genes(adata)


# * Note that: previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA --> do not filter for hv genes here

## Regression & scaling

# convert comma-separate string to list
if ',' in opt["regress_latentvars"]:
    regress_vars = opt["regress_latentvars"].split(',')
else:
    regress_vars = opt["regress_latentvars"]

if opt["regress_latentvars"] == 'none':
    L.info("No regression performed")
else:
    L.info("Starting regression")
    sc.pp.regress_out(adata)

# Clip values exceeding standard deviation 10 according to tutorial.
sc.pp.scale(adata, max_value=10)

# TODO
# run cell cycle scoring
#s_genes = pd.read_csv(opt["sgenes"])[0].tolist()
#g2m_genes = pd.read_csv(opt["g2mgenes"])[0].tolist()
# match the gene ids here!

#sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.write(results_file)


# ########################################################################### #
# ###################### Run PCA & Harmony ################################## #
# ########################################################################### #

## Run PCA
sc.tl.pca(adata, n_comps = opt["nPCs"])

## extract for harmony
if opt["tool"] == 'harmony' :
    L.info("Running harmony")
    data_mat = adata.obsm['X_pca']
    meta_data = adata.obs
    vars_use = [opt["split_var"]]

    ho = hm.run_harmony(data_mat, meta_data, vars_use,
                        sigma = opt["sigma"],
                        plot_convergence = True, max_iter_kmeans=30)

    L.info("Finished harmony")
    adjusted_pcs = pd.DataFrame(ho.Z_corr).T

    adata.obsm['X_harmony']=adjusted_pcs.values

    ## save harmony components
    harmony_out = pd.DataFrame(adata.obsm["X_harmony"],
                               index=adata.obs['barcode'])
    harmony_out.columns = ["harmony_" + str(i) for i in range(1,harmony_out.shape[1]+1)]
    harmony_out.reset_index(inplace=True)

    harmony_out.to_csv(os.path.join(opt["outdir"], "harmony.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

## save anndata object
adata.write(results_file)
L.info("Completed")
