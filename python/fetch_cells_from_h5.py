# This script extracts the data for a list of cells for source h5 files
# from multiple libraries and returns the results as a single anndata.
# Note: the union of the features is returned, with missing feature values
# set to nan values.

import os
import sys
import anndata as ad
import scanpy as sc
import pandas as pd
import logging
import argparse
import numpy as np
import cellhub.tasks.h5 as h5

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
parser.add_argument("--feature_type", default="ALL", type=str,
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

    # the data can be in .h5 [cellranger|cellbender] or .h5ad format.

    h5_path = os.path.join(args.api,
                           "counts","filtered",
                           library_id,
                           "h5","data.h5")
    
    if not os.path.exists(h5_path):
    
        h5_path = os.path.join(args.api,
                           "counts","filtered",
                           library_id,
                           "h5","data.h5ad")

        if not os.path.exists(h5_path):
        
            raise ValueError("h5 file not found on api for: " + library_id)
    
   
    x = h5.read_h5(h5_path)
        
    
           
    # < --------- Extract the requested cells from this library ------------- >
    # 
    common_barcodes = [y for y in barcodes if y in x.obs.index]
     
    L.info("Number of barcodes found:")
    print(len(common_barcodes))

    if len(common_barcodes) == 0:
        raise ValueError("None of the requested barcodes were found in " + h5_path)

    if len(common_barcodes) < len(barcodes):
        L.warn("Not all of the requested barcodes were found in " + h5_path)
        
    L.info("Number of barcodes to extract:")
    print(len(barcodes))
    
    # subset to the set of requested barcodes that are present.          
    x = x[common_barcodes]

    # < ---------------------- Handle the .var frame ------------------------ >

    # use unique identifiers for the index 
    # the index is the ensembl gene name, these are not necessarily unique
    x.var_names_make_unique()
    
    # Raise an error if a sample has no features.
    if len(x.var.index) == 0:
        raise ValueError("No data found in " + h5_path)

    if first:
        # the anndata concat function does not preserve the .var frame (only 
        # index values) so we need to manually keep track of this and then add 
        # the information back after concatenation. 
        var_frame = x.var.copy()
        first = False
        
    else:
        if not all([x for x in x.var.index in var_frame.index]):
            L.info("Adding new .var features from library " + library_id)

            # Update the var_frame to add the new vars.
            var_frame = pd.concat([var_frame.copy(), x.var.copy()],
                                    axis=0).drop_duplicates()          


    # < ---------------------- Handle the .obs frame ------------------------ >

    # It is necessary to construct a unique index for the .var
    x.obs.index = [ y.split("-")[0] + "-" + library_id 
                    for y in x.obs.index.values ]

    anndata_slices[library_id] = x.copy()

L.info("Stitching the slices together")

# Here we perform an outer join to use the union of the
# variables. This is useful when e.g. some samples are missing
# ADT or GEX features. 
#
# Missing features values are populated with
# numpy nan values.
#
anndata = ad.concat(anndata_slices,
                    join="outer",
                    fill_value=np.nan)


L.info("Adding the metadata")
anndata.obs = cell_table.loc[anndata.obs.index,]

# the anndata concate function does not preserve
# the .var frame so we add this back here.
L.info("Adding the var")
anndata.var = var_frame.loc[anndata.var.index]

# drop unnecessary .obs columns
anndata.obs.drop("barcode", axis=1, inplace=True)


if args.feature_type.upper() != "ALL":
    # Subset to the desired feature type.
    if args.feature_type.upper() == "GEX":
        L.info("subsetting to Gene Expression features")
        anndata = anndata[:,anndata.var.feature_types == 'Gene Expression']
    elif args.feature_type.upper() == "ADT":
        L.info("subsetting to ADT features")
        anndata = anndata[:,anndata.var.feature_types == 'Antibody Capture']
    else:
        raise ValueError("Unrecognised feature type")


# save the anndata
anndata.write_h5ad(os.path.join(args.outdir, args.outname))

L.info("complete")
