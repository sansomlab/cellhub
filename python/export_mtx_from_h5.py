import os
import sys
import anndata as ad
import pandas as pd
import logging
import gzip
import argparse
import scipy.io as scio
import cellhub.tasks.h5 as h5

# < ------------------------- Set up the logging --------------------------> #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("extract_cells_from_h5.py")

# < --------------------------- parse arguments ----------------------------> #

parser = argparse.ArgumentParser()
parser.add_argument("--h5", default=None, type=str,
                    help=('the path to a h5 file in cellranger, '
                          'cellbender or anndata format'))
parser.add_argument("--mtx_dir",default=None, type=str,
                    help="output directory for mtx, barcodes and features")

args = parser.parse_args()

L.info("args:")
print(args)

# < ---------------------------- export the mtx ----------------------------> #

if not os.path.exists(args.mtx_dir):
    os.makedirs(args.mtx_dir)

L.info("Reading the h5 file")
        
x = h5.read_h5(args.h5,
               cellbender_analyzed_barcodes_only=False)

# define the output paths
matrix_loc = os.path.join(args.mtx_dir,"matrix.mtx.gz")
barcodes_loc = os.path.join(args.mtx_dir,"barcodes.tsv.gz")
features_loc = os.path.join(args.mtx_dir,"features.tsv.gz")

L.info("Saving the matrix")
with gzip.open(matrix_loc, "w") as matrix_out:
    scio.mmwrite(matrix_out, x.X.T, comment='', 
                 field="integer", precision=None, symmetry=None)
    
L.info("Saving the barcodes")
# columns: barcode #
x.obs.to_csv(barcodes_loc, index=True, columns=[], header=False)

L.info("Saving the features")
# columns: gene_id, gene_name, feature_type

x.var["gene_name"] = x.var.index
feature_cols = ["gene_ids", "gene_name", "feature_types"]

x.var.to_csv(features_loc, index=False, header=False, 
             columns=feature_cols, sep="\t")
  
L.info("complete")
