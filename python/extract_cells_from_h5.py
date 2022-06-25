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
L = logging.getLogger("extract_cells_from_loom.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--cells", default=None, type=str,
                    help='a file containing a mapping of "barcode" to "sequencing_id"')
parser.add_argument("--api", default=None, type=str,
                    help='the path to the api')
#parser.add_argument("--colname", default=None, type=str,
#                    help='column name to use for extraction from samples file')
parser.add_argument("--outdir",default=None, type=str,
                    help="the place to write the final loom")

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

anndata_slices = []

for library_id in libraries:

    # get the barcodes to slice out
    barcodes = cell_table["barcode_id"].values[cell_table["library_id"]==library_id]
    barcodes = [ x.split("_")[0] for x in barcodes ]

    # read h5
    # cellranger.multi/GEX/filtered/${library_id}/h5/sample_feature_bc_matrix.h5
    h5_path = os.path.join(args.apidir,
                           "cellranger.multi","GEX","filtered",
                           library_id,
                           "h5","sample_feature_bc_matrix.h5")

    x = sc.read_10x_h5(h5_path)
    x = x[barcodes]
    x.obs.index = [ y + "_" + library_id for y in x.obs.index.values ]
    x.obs["barcode_id"] = x.obs.index
    x.obs["library_id"] = library_id

    anndata_slices.append(x.copy())


anndata = ad.concat(anndata_slices)
anndata.obs = cell_table.ix[anndata.obs.index,]

anndata.write_h5ad(os.path.join(args.outdir,
"anndata.h5ad"))

L.info("complete")
