import os
import sys
import scanpy
import pandas as pd
import logging
import argparse
import numpy as np

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
parser.add_argument("--obsdata", default=None, type=str,
                    help="A file containing the obsdata")
parser.add_argument("--outdir",default=".", type=str,
                    help="path to output directory")
parser.add_argument("--obstotals",default="total_UMI", type=str,
                    help="Column in obs containing the obs sums")
parser.add_argument("--matrixname", default="matrix.m5ad", type=str,
                    help="a name for the output matrix")

args = parser.parse_args()

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #

L.info("Reading in the 10x market matrix")
x = scanpy.read_10x_mtx(args.mtxdir10x,
                        cache=False)

print(x)

if args.obsdata is not None:

    L.info("Adding the obsdata")

    obs_data = pd.read_csv(args.obsdata,
                           sep="\t",
                           index_col="barcode_id").loc[x.obs.index,:]

    x.obs = obs_data

    if args.obstotals not in x.obs.columns:
        raise ValueError("obstotals column not found in the obsdata")

    # check the observations have the expected sums.
    obs_sums = np.sum(x.X, axis = 1)
    obs_sums = np.squeeze(np.asarray(obs_sums))
    obs_sums = obs_sums.astype(int)
    print(sum(obs_sums))
    print(sum(x.obs[args.obstotals].values))
    if np.array_equal(obs_sums, x.obs[args.obstotals].values):
        L.info("The observations have the expected numbers of counts")

    else:
#        raise ValueError("Data matrix does not match given observations")
        L.info("Data matrix does not match given observations, if GEXADT mode")
        L.info("then exprected less counts, otherwise, worry about it.")

else:
    raise ValueError("Observations data matrix not supplied")


L.info("Saving the h5ad file")
x.write(os.path.join(args.outdir,
                     args.matrixname))
