
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import sys
import os
import logging
import argparse
from datetime import date

from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("tcr_scirpy.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--contig_path", default="", type=str,
                    help="path to the contig annotations")

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

for

ir.io.read_10x_vdj(adata)


ir.tl.chain_qc(adata) #it adds receptor type, subtype, chain pairing
# ########################################################################### #
# ################################### save and exit ######################### #
# ########################################################################### #

L.info("saving the object (with compression")

# with open() as ...:
#   for cell in exp:
#       write(result)

# scirpy_pandas_results_frame.to_csv(arg.outfile, sep="\t", index=False)

L.info("all done")