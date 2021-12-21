#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import os
import re
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

def addMetadata(adata, metadata_infile, id_col):
    metadata_add = pd.read_csv(metadata_infile, sep = "\t",
                                          compression="gzip")
    metadata_add.drop_duplicates(inplace=True)
    new_obs = pd.merge(adata.obs.copy(), metadata_add,
                       on= id_col, how='left')
    new_obs.set_index('barcode_id', inplace=True)
    new_obs = new_obs.loc[adata.obs.index,:]
    adata.obs = new_obs
    return adata

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

L.warning("Running with options ---> %s", opt)

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

# add metadata
adata.obs['barcode_id'] = adata.obs.index.values

if 'metadata_file' in opt.keys():
    adata = addMetadata(adata = adata,
                        metadata_infile = opt['metadata_file'],
                        id_col = opt['metadata_id'])
    L.warning("Added metadata, the following columns are now in the obs: %s",
    str(adata.obs.columns))

# ########################################################################### #
# ################ Run normalization and scaling ############################ #
# ########################################################################### #

## Basic normalization & log-transformation

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

L.warning("Finished normalization and log-transform")

L.warning("The anndata object has %s rows.", adata.shape[0])
L.warning("The anndata object has %s columns.", adata.shape[1])

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


# store the object at this point (with log-normalized but all genes)
if 'metadata_file' in opt.keys():
    adata_full = removeMetadata(adata = adata.copy(),
                                metadata_infile = opt["metadata_file"],
                                id_col = opt["metadata_id"])
else:
    adata_full = adata.copy()

adata_full.write(results_file_logn)

del adata_full

##############################################################################
######################## Identify highly variable genes ######################
##############################################################################

if 'hv_genes' in opt.keys():
    L.warning("Use hv genes from input list")
    gene_list = pd.read_csv(opt["hv_genes"], header=None, sep="\t")
    gene_list.columns = ["gene_name", "gene_id"]
    L.warning("Number of hv genes used: " + len(gene_list) + " genes")
    # make pandas series from hv file
    hvgenes = pd.DataFrame(adata.var["gene_ids"])
    # match hv genes via gene_id
    hvgenes['highly_variable'] = hvgenes['gene_ids'].isin(gene_list['gene_id'])
    hvgenes.columns = ["gene_name", "gene_id", "highly_variable"]
    adata.var['highly_variable'] = hvgenes['highly_variable']
else:
    L.warning("Determine hv genes using scanpy")
    if 'hvg_exclude' in opt.keys():
        L.warning("Remove genes from hv genes as provided here: %s.", opt["hvg_exclude"])
        exclude_list = pd.read_csv(opt["hvg_exclude"], sep="\t")
        L.warning("Number of genes on exclusion list %s genes", len(exclude_list))
        adata = adata[:, ~adata.var['gene_ids'].isin(exclude_list['gene_id'])]
        L.warning("The anndata object has %s rows after excluding genes.", adata.shape[0])
        L.warning("The anndata object has %s columns after excluding genes.", adata.shape[1])
    if opt["merge_hvg"] == "merged":
        sc.pp.highly_variable_genes(adata, n_top_genes=int(opt["ngenes"]))
        hvgenes = adata.var[["highly_variable", "gene_ids"]]
        hvgenes = pd.DataFrame(hvgenes)
    else:
        if ',' in opt["batch_var"]:
            split_vars = opt["batch_var"].split(',')
            hvgenes = adata.var["gene_ids"]
            hvgenes = pd.DataFrame(hvgenes)
            split_vars_cols = []
            for split_var in split_vars:
                split_var_levels = adata.obs[split_var].unique()
                for var_level in split_var_levels:
                    sub_adata = adata[adata.obs[split_var] == var_level].copy()
                    if opt["ngenes"] == "None":
                        sc.pp.highly_variable_genes(sub_adata, batch_key=opt["split_var"])
                        if opt["vote_count"]:
                            maxvotes = sub_adata.obs[opt["split_var"]].nunique()
                            minvotes = round(maxvotes*opt["majority_fraction"])
                            hvgenes.loc[:, 'highly_variable' + '_' + split_var + '_' + var_level] = np.where(sub_adata.var["highly_variable_nbatches"] >= minvotes, True, False)
                            hvgenes.loc[:, 'highly_variable_nbatches' + '_' + split_var + '_' + var_level] = sub_adata.var['highly_variable_nbatches'].copy()
                    else:
                        sc.pp.highly_variable_genes(sub_adata, batch_key=opt["split_var"], n_top_genes=int(opt["ngenes"]))
                        hvgenes.loc[:, 'highly_variable_nbatches' + '_' + split_var + '_' + var_level] = sub_adata.var['highly_variable_nbatches'].copy()
                        hvgenes.loc[:, 'highly_variable' + '_' + split_var + '_' + var_level] = sub_adata.var['highly_variable'].copy()
                    split_vars_cols.append('highly_variable' + '_' + split_var + '_' + var_level)
            hvgenes['highly_variable'] = hvgenes.loc[:, split_vars_cols].transpose().any()
            sc.pp.highly_variable_genes(adata, batch_key=opt["split_var"])
            adata.var.loc[:, 'highly_variable'] = hvgenes['highly_variable'].copy()
        else:
            if opt["ngenes"] == "None":
                sc.pp.highly_variable_genes(adata, batch_key=opt["split_var"])
            else:
                sc.pp.highly_variable_genes(adata, batch_key=opt["split_var"], n_top_genes=int(opt["ngenes"]))
            # extract hv genes and store them
            hvgenes = adata.var[["highly_variable", "gene_ids", "highly_variable_nbatches"]]
            hvgenes = pd.DataFrame(hvgenes)
    
    L.warning("hvg: " + str(hvgenes))
    L.warning("hvg colummns: " + str(hvgenes.columns))
    L.warning(str(adata.var.index))
    L.warning(str(adata.var.columns))
 
    hvgenes.loc[:, "gene_name"] = hvgenes.index.copy()
    hvgenes = hvgenes.loc[hvgenes['highly_variable'] == True]
    hvgenes.reset_index(inplace=True)
    #hvgenes.columns = ["gene_name", "highly_variable", "gene_id"]
    #del hvgenes["highly_variable"]
    hvgenes.to_csv(os.path.join(opt["outdir"], "hv_genes.tsv.gz"),
                       sep="\t", index=False, compression="gzip")
    
sc.pl.highly_variable_genes(adata, save = "_hvg.pdf", show=False)

# make log-norm layer for later use
adata.layers['log1p'] = adata.X.copy()

# add hvg back to full object
full_adata = anndata.read_h5ad(results_file_logn, backed="r+")
full_adata.var['highly_variable'] = full_adata.var['gene_ids'].isin(hvgenes['gene_ids'])
full_adata.write()


if opt["regress_cellcycle"] == 'none':
    ## subset to hv genes for all downstream analyses
    L.warning("Subset object to hv genes only")
    adata = adata[:, adata.var.highly_variable]
else:
    L.warning("Cell cycle regression is used, therefore hv + cell cycle genes are kept. Please note that this will take longer than only using hv genes.")
    keep_genes = list(set(adata.var.loc[adata.var.highly_variable].index.to_list() + cell_cycle_genes))
    adata = adata[:, keep_genes]

L.warning("The anndata object has %s rows after subsetting.", adata.shape[0])
L.warning("The anndata object has %s columns after subsetting.", adata.shape[1])

##################################################################################
################################## Regression & scaling ##########################
##################################################################################
# if more than one latent variable, convert comma-separate string to list
if ',' in opt["regress_latentvars"]:
    regress_vars = opt["regress_latentvars"].split(',')
else:
    if opt["regress_latentvars"] != 'none':
        regress_vars = [opt["regress_latentvars"]]

if opt["regress_latentvars"] == 'none':
    L.warning("No regression performed")
    if opt["regress_cellcycle"] != 'none':
        L.warning("Cell cycle scoring performed. It is set to " + str(opt["regress_cellcycle"]))
        if opt["regress_cellcycle"] == "all":
            regress_vars = ['S_score', 'G2M_score']
        elif opt["regress_cellcycle"] == "difference":
            regress_vars = ['CC.Difference']
        else:
            raise Exception('Cell cycle option not supported. Use all or difference.')
        L.warning("Full list of variables for regression is: " + ",".join(regress_vars))
        L.warning("Starting regression")
        sc.pp.regress_out(adata, regress_vars)
else:
    if opt["regress_cellcycle"] != 'none':
        L.warning("Cell cycle scoring performed. It is set to " + str(opt["regress_cellcycle"]))
        if opt["regress_cellcycle"] == "all":
            regress_vars = regress_vars + ['S_score', 'G2M_score']
        elif opt["regress_cellcycle"] == "difference":
            regress_vars = regress_vars + ['CC.Difference']
        else:
            raise Exception('Cell cycle option not supported. Use all or difference.')
    L.warning("Full list of variables for regression is: " + ",".join(regress_vars))
    L.warning("Starting regression")
    sc.pp.regress_out(adata, regress_vars)

# Clip values exceeding standard deviation 10 according to tutorial.
sc.pp.scale(adata, max_value=10)

if opt["regress_cellcycle"] != 'none':
    L.warning("Check cell cycle effects post-regression")
    adata_cc_genes = adata[:, adata.var.index.isin(cell_cycle_genes)].copy()
    L.warning("Check cell cycle effects adata %s columns. ", adata_cc_genes.shape[1])
    sc.tl.pca(adata_cc_genes)
    sc.pl.pca_scatter(adata_cc_genes, color='phase', show=False, save = "_cc_phase_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='G2M_score', show=False, save = "_cc_G2Mscore_postCorr.pdf")
    sc.pl.pca_scatter(adata_cc_genes, color='S_score', show=False, save = "_cc_Sscore_postCorr.pdf")
    L.warning("Subset object to only hv genes")
    L.warning("hvg: " + str(hvgenes))
    L.warning("hvg colummns: " + str(hvgenes.columns))
    L.warning(str(adata.var.index))
    L.warning(str(adata.var.columns))
    adata = adata[:, adata.var.index.isin(hvgenes["gene_name"])]
    L.warning("The anndata object has %s columns after subsetting for hv genes.", adata.shape[1])
else:
    L.warning("Cell cycle scoring set to none.")

adata.write(results_file)


# ########################################################################### #
# ###################### Run PCA & Harmony ################################## #
# ########################################################################### #

## Run PCA
sc.tl.pca(adata, n_comps = opt["totalPCs"])
sc.pl.pca_variance_ratio(adata, save = "_pca_stdev.pdf",
                         show=False, n_pcs = opt["totalPCs"])
sc.pl.pca_variance_ratio(adata, save = "_pca_stdev_log.pdf",log=True,
                         show=False, n_pcs = opt["totalPCs"])
sc.pl.pca_loadings(adata, show=False,save="_pca_top_loadings.pdf")

# write out the variance ratios
pca_var = pd.DataFrame(adata.uns['pca']['variance_ratio'])
pca_var.columns = ['pca_variance_ratio']
pca_var.to_csv(os.path.join(opt["outdir"], "pca_variance_ratios.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

# write out PCA loadings
df_loadings = pd.DataFrame(adata.varm['PCs'], index=adata.var_names)
df_loadings.index.name = 'gene_name'
df_loadings.reset_index(inplace=True)
df_loadings.to_csv(os.path.join(opt["outdir"], "pca_loadings.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

## extract for harmony
if opt["tool"] == 'harmony' :
    L.warning("Running harmony")
    data_mat = adata.obsm['X_pca'][:, 0:opt["nPCs"]].copy()
    if ',' in opt["split_var"]:
        vars_use = opt["split_var"].split(',')
    else:
        vars_use = [opt["split_var"]]

    meta_data = adata.obs[vars_use]

    L.warning("Using following variables for harmony integration: " + ",".join(vars_use))

    ho = hm.run_harmony(data_mat, meta_data, vars_use,
                        sigma = opt["sigma"],
                        #plot_convergence = True, # deprecated
                        max_iter_kmeans=30,
                        max_iter_harmony = opt["max_iter_harmony"])


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

    ## save harmony objective values (convergence visualization)
    harmony_obj = pd.DataFrame(ho.objective_harmony.copy())
    harmony_obj.index.name = "iteration"
    harmony_obj.columns = ['objective']
    harmony_obj.reset_index(inplace=True)
    harmony_obj.to_csv(os.path.join(opt["outdir"], "harmony_objective.tsv.gz"),
                       sep="\t", index=False, compression="gzip")

elif opt["tool"] == "bbknn":
    L.warning("Running bbknn, no dim reduction will be stored")
    L.warning("Using %s PCA components for bbknn", opt["nPCs"])
    bbknn.bbknn(adata, batch_key=opt["split_var"], n_pcs = opt["nPCs"])

elif opt["tool"] == "scanorama":
    L.warning("Splitting anndata object for scanorama")
    list_ids = (set(adata.obs[opt["split_var"]].tolist()))
    all_anndata = []

    for p in list_ids:
        all_anndata = all_anndata + [adata[adata.obs[opt["split_var"]] == p]]

    # The batch size used in the alignment vector computation. Useful when
    # correcting very large (>100k samples) data sets. Set to large value
    # that runs within available memory.
    # hvg= option could be used but object already subset to hvg
    integrated = scanorama.integrate_scanpy(all_anndata, dimred=opt["nPCs"], batch_size=5000)

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
    pca_out = pd.DataFrame(adata.obsm['X_pca'][:, 0:opt["nPCs"]].copy(),
                           index = adata.obs.index.copy())
    pca_out.index.name = 'barcode'
    pca_out.columns = ["PC_" + str(i) for i in range(1,pca_out.shape[1]+1)]
    pca_out.reset_index(inplace=True)
    pca_out.to_csv(os.path.join(opt["outdir"], "pca.tsv.gz"),
                   sep="\t", index=False, compression="gzip")

## save anndata object
adata.write(results_file)
L.warning("Completed")
