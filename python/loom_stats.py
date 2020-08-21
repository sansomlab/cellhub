import os
import sys
import loompy
import pandas as pd
import logging
import argparse
import numpy as np
from scipy.stats import spearmanr

# compute some basic stats on the cells present in a loom file for
# sanity checking.

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("loom_stats.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--loom", default=None, type=str,
                    help='a loom file')
parser.add_argument("--reference", default=None, type=str,
                     help='a table containing pre-computed reference statistics')
parser.add_argument("--compare", default=None, type=str,
                     help='a comma separated list of the statistics to comapare with the reference')
parser.add_argument("--outdir",default=None, type=str,
                    help="the place to write the stats loom")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################## read in the matrix ######################### #
# ########################################################################### #

L.info("computing stats")

with loompy.connect(args.loom) as ds:

    # velocyto.py mangles the barcodes!
    # from AAAAAA-1
    # to   sample:AAAAAx

    total_counts = ds.map([np.sum], axis=1)[0]

    barcodes = ds.ca['CellID']

    # print(ds.ca['CellID'][0:10])

    # this is not good.. input files are being modified!


L.info("saving the computed stats")

results = pd.DataFrame({"total_UMI": total_counts,
                        "barcode_id": barcodes})

results.to_csv(os.path.join(args.outdir,
                            "loom.stats.tsv.gz"),
               sep="\t",
               index=False)


L.info("comparing with reference")

reference = pd.read_csv(args.reference, sep="\t")

stats = [x.strip() for x in args.compare.split(",")]

for stat in stats:

    if stat == "barcode_id":
        # check (1) memebership and (2) order

        x = reference.barcode_id
        y = results.barcode_id

        xy = list(set(x) & set(y))

        L.info("no. barcodes in reference: " + str(len(x)))
        L.info("no. barcodes in loom: " + str(len(y)))
        L.info("no. barcodes in both: " + str(len(xy)))


    if stat == "total_UMI":

        x = reference.total_UMI
        y = results.total_UMI

        L.info("total counts in reference:" + str(np.sum(x)))
        L.info("total counts in loom: " + str(np.sum(y)))


        reference.index = reference.barcode_id
        x.ordered = reference.loc[results.barcode_id,"total_UMI"]

        L.info("correlation in total umi with unsorted barcode order: " + str(spearmanr(x, y)))
        L.info("correlation in total umi with matching barcode order: " + str(spearmanr(x.ordered, y)))


L.info("complete")
