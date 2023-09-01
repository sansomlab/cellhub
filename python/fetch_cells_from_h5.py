import os
import sys
import anndata as ad
import scanpy as sc
import pandas as pd
import logging
import argparse
import numpy as np
import cellhub.tasks.cellbender as cb

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("extract_cells_from_h5.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--cells", default=None, type=str,
                    help='a file containing a mapping of "barcode" to "sequencing_id"')
parser.add_argument("--api", default=None, type=str,
                    help='the path to the api')
parser.add_argument("--feature_type", default=None, type=str,
                    help='for "Gene Expression" use "GEX", for "Antibody Capture" use "ADT"')
parser.add_argument("--outdir",default=None, type=str,
                    help="the output directory")
parser.add_argument("--outname",default=None, type=str,
                    help="the output file name")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #

L.info("Reading cell table")

cell_table = pd.read_csv(args.cells, sep="\t")

good_cols = [x for x in cell_table.columns if
             not x.startswith("barcode:") and
             not x.startswith("library_id:") and
             not x.startswith("filename")]

cell_table = cell_table[good_cols]

cell_table.index = ["-".join(y) for y in 
                    zip([x.split("-")[0] for x in cell_table["barcode"].values],
                        cell_table["library_id"].values)]  

libraries = set(cell_table["library_id"].values)

anndata_slices = {}

first = True
for library_id in libraries:

    L.info("Slicing " + library_id)

    # get the barcodes to slice out
    barcodes = cell_table["barcode"].values[cell_table["library_id"]==library_id]

    h5_path = os.path.join(args.api,
                           "counts","filtered",
                           library_id,
                           "h5","data.h5")

    if not os.path.exists(h5_path):
        raise ValueError("h5 path does not exist: " + h5_path)

    try:
        
        x = sc.read_10x_h5(h5_path)

    except:
        # we need to use a custom loader for cellbender
        L.info("Reading the h5 file with scanpy read_10x_h5 failed... "
               "attempting to read with the custom cellbender loader")
    
        h5_path = os.path.join(args.api,
                           "counts","filtered",
                           library_id,
                           "h5","data.h5")
        
        x = cb.anndata_from_h5(h5_path)
        
        if "feature_type" in x.var.columns:
        
            x.var.columns = [y.replace("feature_type","feature_types") for y in x.var.columns]
        
        if "gene_id" in x.var.columns:

            x.var.columns = [y.replace("gene_id","gene_ids") for y in x.var.columns]
        

    
    # Subset to the desired feature type.
    if args.feature_type == "GEX":
        L.info("subsetting to Gene Expression features")
        x = x[:,x.var.feature_types == 'Gene Expression']
    elif args.feature_type == "ADT":
        L.info("subsetting to Gene Expression features")
        x = x[:,x.var.feature_types == 'Antibody Capture']
    else:
        if x.args.feature_type is not None:
            raise ValueError("Unrecognised feature type")

    # use unique identifiers for the index 
    # the index is the ensembl gene name, these are not necessarily unique
    x.var_names_make_unique()

    if first:
        var_frame = x.var.copy()
        var_index = x.var.index.copy()
        first = False
    else:
        if not np.array_equal(x.var.index, var_index):
            raise ValueError("var index mismatch")

    L.info("Number of barcodes to extract:")
    print(len(barcodes))
    
    common_barcodes = [y for y in barcodes if y in x.obs.index]
     
    L.info("Number of barcodes found:")
    print(len(common_barcodes))

    if len(common_barcodes) < len(barcodes):
        if args.source == "cellranger":
            raise ValueError("Not all barcodes were found")
        else:
            L.warn("Not all barcodes were found")

    x = x[common_barcodes]

    x.obs.index = [ y.split("-")[0] + "-" + library_id 
                    for y in x.obs.index.values ]

    anndata_slices[library_id] = x.copy()

L.info("Stitching the slices together")
anndata = ad.concat(anndata_slices)


L.info("Adding the metadata")
anndata.obs = cell_table.loc[anndata.obs.index,]

# not sure why this is necessary...
L.info("Adding the var")
anndata.var = var_frame.loc[anndata.var.index]

# drop unnecessary .obs columns
anndata.obs.drop("barcode", axis=1, inplace=True)
    
# save the anndata
anndata.write_h5ad(os.path.join(args.outdir, args.outname))

L.info("complete")
