#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import os
import scanpy as sc
import anndata
import bbknn
import scanorama
import numpy as np
import scanpy.external as sce
import harmonypy as hm
import pandas as pd
import matplotlib.pyplot as plt
import logging
import yaml
import warnings
warnings.filterwarnings('ignore')

# ########################################################################### #
# ############### Set up the log and figure folder ########################## #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.WARNING)
L = logging.getLogger("run_integration")

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

# write two results files 
results_file_logn = os.path.join(opt["outdir"], "lognormalized_anndata.h5ad")
results_file = os.path.join(opt["outdir"], "normalized_integrated_anndata.h5ad")

L.info("Running with options ---> %s", opt)

L.warning("Writing output to file %s", results_file)

# checks
if opt["tool"] == "scanorama" and ',' in opt["split_var"]:
    raise Exception("Scanorama can only take one variable for integration.")

if opt["tool"] == "bbknn" and ',' in opt["split_var"]:
    raise Exception("bbknn can only take one variable for integration.")


# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #

adata = anndata.read(os.path.join(opt["matrixdir"], "matrix.h5ad"))
# create a raw counts layer
adata.layers['counts'] = adata.X.copy()

L.warning("Loaded anndata object")

# write out gene annotation
anno_genes = pd.DataFrame(adata.var["gene_ids"].copy())
anno_genes.reset_index(inplace=True)
anno_genes.columns = ["gene_name", "gene_id"]
anno_genes.to_csv(os.path.join(opt["outdir"], "annotation_genes.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

# drop column named 'barcode' if it exists
if 'barcode' in adata.obs.columns:
    L.warning("Removing column called barcode from .obs, use index")
    del adata.obs['barcode']

# ########################################################################### #
# ################ Run normalization and scaling ############################ #
# ########################################################################### #

## Basic normalization & log-transformation

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

L.warning("Finished normalization and log-transform")

if 'sgenes' in opt.keys() and 'g2mgenes' in opt.keys():
    L.warning("Making separate object for cell cycle scoring")

    adata_cc_genes = adata.copy()
    sc.pp.scale(adata_cc_genes)
    # cell cycle scoring similar to Seurat method
    # see here: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
    L.warning("Running cell cycle scoring")

    s_genes = pd.read_csv(opt["sgenes"], header=None)[0].tolist()             
    g2m_genes = pd.read_csv(opt["g2mgenes"], header=None)[0].tolist()
    
    # check that genes are in adata
    s_genes = [x for x in s_genes if x in adata_cc_genes.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata_cc_genes.var_names]
    L.warning("Number of G2M genes: " + str(len(g2m_genes)))
    L.warning("Number of S genes: " + str(len(s_genes)))

    sc.tl.score_genes_cell_cycle(adata_cc_genes, s_genes=s_genes, g2m_genes=g2m_genes)
    # add score for difference between phases
    adata_cc_genes.obs['CC.Difference'] = adata_cc_genes.obs['S_score'] - adata_cc_genes.obs['G2M_score']

    # assign cell cycle scores back to original object
    adata.obs['S_score'] = adata_cc_genes.obs['S_score']
    adata.obs['G2M_score'] = adata_cc_genes.obs['G2M_score']
    adata.obs['phase'] = adata_cc_genes.obs['phase']
    adata.obs['CC.Difference'] = adata_cc_genes.obs['CC.Difference']

    ## make new object to run PCA ONLY on cc genes (otherwise all genes are used by Scanpy)
    L.warning("Running PCA based on only cell cycle genes")
    cell_cycle_genes = s_genes + g2m_genes
    adata_cc_subset = adata_cc_genes[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_subset)
    sc.pl.pca_scatter(adata_cc_subset, color='phase', show=False, save = "_cc_phase.pdf")
    sc.pl.pca_scatter(adata_cc_subset, color='G2M_score', show=False, save = "_cc_G2Mscore.pdf")
    sc.pl.pca_scatter(adata_cc_subset, color='S_score', show=False, save = "_cc_Sscore.pdf")


## Identify highly variable genes

# If n_top_genes is used then the flavour = 'cellranger'
if 'hv_genes' in opt.keys():
    L.warning("Use hv genes from input list")
    gene_list = pd.read_csv(opt["hv_genes"])
    gene_list = pd.read_csv(opt["hv_genes"], header=None, sep="\t")
    gene_list.columns = ["gene_name", "gene_id"]
    L.warning("Number of hv genes used: " + len(gene_list) + " genes")
    # make pandas series from hv file
    hvgenes = pd.DataFrame(adata.var["gene_ids"])
    # match hv genes via gene_id 
    hvgenes['highly_variable'] = hvgenes['gene_ids'].isin(gene_list['gene_id'])
    adata.var['highly_variable'] = hvgenes['highly_variable'] 
else:
    L.warning("Determine hv genes using scanpy")
    sc.pp.highly_variable_genes(adata, n_top_genes=opt["ngenes"])
    # extract hv genes and store them
    hvgenes = adata.var[["highly_variable","gene_ids"]]
    hvgenes = pd.DataFrame(hvgenes)
    hvgenes = hvgenes.loc[hvgenes['highly_variable'] == True] 
    hvgenes.reset_index(inplace=True)
    hvgenes.columns = ["gene_name", "highly_variable", "gene_id"]
    del hvgenes["highly_variable"]
    hvgenes.to_csv(os.path.join(opt["outdir"], "hv_genes.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

sc.pl.highly_variable_genes(adata, save = "_hvg.pdf", show=False)


# store the objects at this point (with log-normalized and hvg)
adata.write(results_file_logn)
adata.write(results_file)
# make log-norm layer for later use
adata.layers['log1p'] = adata.X.copy()
adata.write(results_file)

## subset to hv genes for all downstream analyses
L.warning("Subset object to hv genes only")
adata = adata[:, adata.var.highly_variable]

## Regression & scaling

# convert comma-separate string to list
if ',' in opt["regress_latentvars"]:
    regress_vars = opt["regress_latentvars"].split(',')
else:
    regress_vars = [opt["regress_latentvars"]]

if opt["regress_latentvars"] == 'none':
    L.warning("No regression performed")
else:
    L.warning("Starting regression")
    sc.pp.regress_out(adata, regress_vars)

# Clip values exceeding standard deviation 10 according to tutorial.
sc.pp.scale(adata, max_value=10)

if opt["regress_cellcycle"] != 'none':
    L.warning("Cell cycle scoring performed. It is set to " + str(opt["regress_cellcycle"]))
    if opt["regress_latentvars"] == 'none':
        regress_vars = []
    elif ',' in opt["regress_latentvars"]:
        regress_vars = opt["regress_latentvars"].split(',')
    else:
        regress_vars = [opt["regress_latentvars"]]
    if opt["regress_cellcycle"] == "all":
        regress_vars = regress_vars + ['S_score', 'G2M_score']
    elif opt["regress_cellcycle"] == "difference":
        regress_vars = regress_vars + ['CC.Difference']
    else:
        raise Exception('Cell cycle option not supported. Use all or difference.')
    L.warning("Full list of variables for regression is: " + ",".join(regress_vars))
    sc.pp.regress_out(adata, regress_vars)
    sc.pp.scale(adata, max_value=10)

    adata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)
    sc.pl.pca_scatter(adata_cc_genes, color='phase', show=False, save = "_cc_phase_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='G2M_score', show=False, save = "_cc_G2Mscore_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='S_score', show=False, save = "_cc_Sscore_postCorr.pdf")
else:
    L.warning("Cell cycle scoring set to none.")

adata.write(results_file)


# ########################################################################### #
# ###################### Run PCA & Harmony ################################## #
# ########################################################################### #

## Run PCA
sc.tl.pca(adata, n_comps = opt["nPCs"])
sc.pl.pca_variance_ratio(adata, save = "_pca_stdev.pdf", 
                         show=False, n_pcs = opt["nPCs"])
sc.pl.pca_variance_ratio(adata, save = "_pca_stdev_log.pdf",log=True, 
                         show=False, n_pcs = opt["nPCs"])

## extract for harmony
if opt["tool"] == 'harmony' :
    L.warning("Running harmony")
    data_mat = adata.obsm['X_pca']
    if ',' in opt["split_var"]:
        vars_use = opt["split_var"].split(',')
    else:
        vars_use = [opt["split_var"]]

    meta_data = adata.obs[vars_use]

    L.warning("Using following variables for harmony integration: " + ",".join(vars_use))

    ho = hm.run_harmony(data_mat, meta_data, vars_use,
                        sigma = opt["sigma"],
                        plot_convergence = True, max_iter_kmeans=30)

    L.warning("Finished harmony")
    adjusted_pcs = pd.DataFrame(ho.Z_corr).T

    adata.obsm['X_harmony']=adjusted_pcs.values
 
    ## save harmony components
    harmony_out = pd.DataFrame(adata.obsm["X_harmony"].copy(),
                               index=adata.obs.index.copy())
    harmony_out.index.name = 'barcode'
    harmony_out.columns = ["harmony_" + str(i) for i in range(1,harmony_out.shape[1]+1)]
    harmony_out.reset_index(inplace=True)

    harmony_out.to_csv(os.path.join(opt["outdir"], "harmony.tsv.gz"),
                       sep="\t", index=False, compression="gzip")
elif opt["tool"] == "bbknn":
    L.warning("Running bbknn, no dim reduction will be stored")
    bbknn.bbknn(adata, batch_key=opt["split_var"], n_pcs = opt["nPCs"])

elif opt["tool"] == "scanorama":
    L.warning("Splitting anndata object for scanorama")
    list_ids = (set(adata.obs[opt["split_var"]].tolist()))
    all_anndata = []

    for p in list_ids:
        all_anndata = all_anndata + [adata[adata.obs[opt["split_var"]] == p]]
    
    #     batch_size: `int`, optional (default: `5000`)
    #         The batch size used in the alignment vector computation. Useful when
    #         correcting very large (>100k samples) data sets. Set to large value
    #         that runs within available memory.
    # batch_size: `int`, optional (default: `5000`)
    #     The batch size used in the alignment vector computation. Useful when
    #     correcting very large (>100k samples) data sets. Set to large value
    #     that runs within available memory.
    ## hvg= option could be used but object already subset to hvg
    integrated = scanorama.integrate_scanpy(all_anndata, dimred=50, batch_size=5000)

    embedding_scanorama = np.concatenate(integrated, axis=0)
    adata.obsm["X_scanorama_embedding"] = embedding_scanorama
    L.warning("Finished scanorama and written embeddings into anndata")
    ## save scanorama components
    scanorama_out = pd.DataFrame(adata.obsm["X_scanorama_embedding"].copy(),
                                 index=adata.obs.index.copy())
    scanorama_out.index.name = 'barcode'
    scanorama_out.columns = ["scanorama_" + str(i) for i in range(1,scanorama_out.shape[1]+1)]
    scanorama_out.reset_index(inplace=True)

    scanorama_out.to_csv(os.path.join(opt["outdir"], "scanorama.tsv.gz"),
                       sep="\t", index=False, compression="gzip")
    L.warning("Finished writing scanorama embeddings to file")
else:
    L.warning("Write out PCA components")
    pca_out = pd.DataFrame(adata.obsm['X_pca'].copy(),
                           index = adata.obs.index.copy())
    pca_out.index.name = 'barcode'
    pca_out.columns = ["PC_" + str(i) for i in range(1,pca_out.shape[1]+1)]
    pca_out.reset_index(inplace=True)
    pca_out.to_csv(os.path.join(opt["outdir"], "pca.tsv.gz"),
                   sep="\t", index=False, compression="gzip")

## save anndata object
adata.write(results_file)
L.warning("Completed")
