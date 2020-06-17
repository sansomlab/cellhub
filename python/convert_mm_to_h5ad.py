import os
import sys
import scanpy
import pandas as pd
import logging
import argparse


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("convert_mtx_to_h5ad.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--mtxdir10x", default=".", type=str,
                    help="the directory containing the 10x matrix")
parser.add_argument("--metadata", default=None, type=str,
                    help="A file containing the metadata")
parser.add_argument("--qcdata", default=None, type=str,
                    help="A file containing the qc info")
parser.add_argument("--outdir",default=".", type=str,
                    help="path to output directory")
parser.add_argument("--matrixname", default="matrix.m5ad", type=str,
                    help="a name for the output matrix")

args = parser.parse_args()

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #

L.info("Reading in the 10x market matrix")
x = scanpy.read_10x_mtx(args.mtxdir10x,
                        cache=True)

if args.metadata is not None:

    L.info("Adding the metadata")

    x.uns["metadata"] = pd.read_csv(args.metadata, 
                                    sep="\t",
                                    index_col="barcode").loc[x.obs.index,:]


if args.qcdata is not None:

    L.info("Adding the qc data")

    x.uns["qcdata"] = pd.read_csv(args.qcdata, 
                                    sep="\t",
                                    index_col="barcode").loc[x.obs.index,:]


L.info("Saving the h5ad file")
x.write(os.path.join(args.outdir,
                     args.matrixname))
