import numpy as np
import pandas as pd
import scipy as sc
### import scirpy ?!?
import sys
import os
import logging
import argparse
from datetime import date


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("tcr_scirpy.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--cellranger_vdj_ouput", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")

parser.add_argument("--param1", default="umap.tsv.gz", type=str,
                    help="a comma seperated list of umaps to add")

parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################### Script.. ################################# #
# ########################################################################### #


L.info("Running scirpy")

# code to run scirpy goes here...

# ########################################################################### #
# ################################### save and exit ######################### #
# ########################################################################### #

L.info("saving the object (with compression")

# with open() as ...:
#   for cell in exp:
#       write(result)

# scirpy_pandas_results_frame.to_csv(arg.outfile, sep="\t", index=False)

L.info("all done")
