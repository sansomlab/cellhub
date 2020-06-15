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
sc.settings.figdir = os.path.join(opt["outdir"], "figures.dir")
sc.settings.set_figure_params(dpi=300, dpi_save=300)

# write folder                                                                                    
results_file = os.path.join(opt["outdir"], "normalized_integrated_anndata.h5ad")

L.info("Running with options ---> %s", opt)

L.info("Writing output to file %s", results_file)

# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #


adata = sc.read_10x_mtx(opt["matrixdir"],
                        var_names='gene_symbols')
# could use gene_symbols here but then also required:
adata.var_names_make_unique()

# ## add metadata for cells
metadata = pd.read_csv(os.path.join(opt["matrixdir"], "metadata.tsv.gz"), sep="\t")
adata.obs = metadata

L.info("Made anndata object")

# write out gene annotation
anno_genes = pd.DataFrame(adata.var["gene_ids"])
anno_genes.reset_index(inplace=True)
anno_genes.columns = ["gene_name", "gene_id"]
anno_genes.to_csv(os.path.join(opt["outdir"], "annotation_genes.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

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
    L.info("Use hv genes from input list")
    gene_list = pd.read_csv(opt["hv_genes"])
    gene_list = pd.read_csv(opt["hv_genes"], header=None, sep="\t")
    gene_list.columns = ["gene_name", "gene_id"]
    L.info("Number of hv genes used: " + len(gene_list) + " genes")
    # make pandas series from hv file
    hvgenes = pd.DataFrame(adata.var["gene_ids"])
    # match hv genes via gene_id 
    hvgenes['highly_variable'] = hvgenes['gene_ids'].isin(gene_list['gene_id'])
    adata.var['highly_variable'] = hvgenes['highly_variable'] 
else:
    L.info("Determine hv genes using scanpy")
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

# * Note that: previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA --> do not filter for hv genes here

## Regression & scaling

# convert comma-separate string to list
if ',' in opt["regress_latentvars"]:
    regress_vars = opt["regress_latentvars"].split(',')
else:
    regress_vars = [opt["regress_latentvars"]]

if opt["regress_latentvars"] == 'none':
    L.info("No regression performed")
else:
    L.info("Starting regression")
    sc.pp.regress_out(adata, regress_vars)

# Clip values exceeding standard deviation 10 according to tutorial.
sc.pp.scale(adata, max_value=10)

if 'sgenes' in opt.keys() and 'g2mgenes' in opt.keys():
    # cell cycle scoring similar to Seurat method
    # see here: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
    L.info("Running cell cycle scoring")

    s_genes = pd.read_csv(opt["sgenes"], header=None)[0].tolist()             
    g2m_genes = pd.read_csv(opt["g2mgenes"], header=None)[0].tolist()
    
    # check that genes are in adata
    s_genes = [x for x in s_genes if x in adata.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata.var_names]
    L.info("Number of G2M genes: " + str(len(g2m_genes)))
    L.info("Number of S genes: " + str(len(s_genes)))

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # add score for difference between phases
    adata.obs['CC.Difference'] = adata.obs['S_score'] - adata.obs['G2M_score']

    ## make new object to run PCA ONLY on cc genes (otherwise all genes are used by Scanpy)
    L.info("Running PCA based on only cell cycle genes")
    cell_cycle_genes = s_genes + g2m_genes
    adata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)
    sc.pl.pca_scatter(adata_cc_genes, color='phase', show=False, save = "_cc_phase.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='G2M_score', show=False, save = "_cc_G2Mscore.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='S_score', show=False, save = "_cc_Sscore.pdf")


if opt["regress_cellcycle"] != 'none':
    L.info("Cell cycle scoring performed. It is set to " + str(opt["regress_cellcycle"]))
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
    L.info("Full list of variables for regression is: " + ",".join(regress_vars))
    sc.pp.regress_out(adata, regress_vars)
    sc.pp.scale(adata, max_value=10)

    adata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)
    sc.pl.pca_scatter(adata_cc_genes, color='phase', show=False, save = "_cc_phase_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='G2M_score', show=False, save = "_cc_G2Mscore_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='S_score', show=False, save = "_cc_Sscore_postCorr.pdf")
else:
    L.info("Cell cycle scoring set to none.")


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
    if ',' in opt["split_var"]:
        vars_use = opt["split_var"].split(',')
    else:
        vars_use = [opt["split_var"]]

    L.info("Using following variables for harmony integration: " + ",".join(vars_use))

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
else:
    L.info("Write out PCA components")
    pca_out = pd.DataFrame(adata.obsm['X_pca'])
    pca_out.columns = ["PC_" + str(i) for i in range(1,pca_out.shape[1]+1)]
    pca_out['barcode'] = adata.obs['barcode'].tolist()
    pca_out.to_csv(os.path.join(opt["outdir"], "pca.tsv.gz"),
                   sep="\t", index=False, compression="gzip")

## save anndata object
adata.write(results_file)
L.info("Completed")
