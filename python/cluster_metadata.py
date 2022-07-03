import os
import re
import argparse
import anndata as ad
import pandas as pd
import logging
import sys



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output metadata file ")
args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# ########################################################################### #
# #########################  Make the anndata object ######################### #
# ########################################################################### #

adata = ad.read_h5ad(args.source_anndata)

# compute clusters
adata.obs.to_csv(args.outfile, sep="\t", index=False)                       

L.info("Complete")
