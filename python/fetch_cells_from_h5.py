import os
import sys
import anndata as ad
import scanpy as sc
import pandas as pd
import logging
import argparse
import numpy as np

#    statement = '''python %(cellhub_dir)s/python/extract_cells_from_h5.py"
#                       --cells=%(cell_manifest)s
#                       --samples=%(sample_table)s
#                       --outdir=%(outdir)s
#                '''

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
cell_table.index = cell_table["barcode_id"].values

libraries = set(cell_table["library_id"].values)

anndata_slices = {}

first = True
for library_id in libraries:

    L.info("Slicing " + library_id)

    # get the barcodes to slice out
    barcodes = cell_table["barcode_id"].values[cell_table["library_id"]==library_id]
    barcodes = [ x.split("-")[0] + "-1" for x in barcodes ]

    h5_path = os.path.join(args.api,
                           "cellranger.multi","counts","filtered",
                           library_id,
                           "h5","sample_filtered_feature_bc_matrix.h5")

    x = sc.read_10x_h5(h5_path)
    
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

    x = x[barcodes]
    x.obs.index = [ y.split("-")[0] + "-" + library_id 
                    for y in x.obs.index.values ]
    x.obs["barcode_id"] = x.obs.index
    x.obs["library_id"] = library_id

    anndata_slices[library_id] = x.copy()

L.info("Stitching the slices together")
anndata = ad.concat(anndata_slices)


L.info("Adding the metadata")
anndata.obs = cell_table.loc[anndata.obs.index,]

# not sure why this is necessary...
L.info("Adding the var")
anndata.var = var_frame.loc[anndata.var.index]

anndata.write_h5ad(os.path.join(args.outdir, args.outname))

L.info("complete")
