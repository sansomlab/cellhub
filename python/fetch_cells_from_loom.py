import os
import sys
import loompy
import pandas as pd
import logging
import argparse
import numpy as np

#    statement = '''python %(cellhub_dir)s/python/extract_cells_from_loom.py"
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
parser.add_argument("--samples", default=None, type=str,
                    help='a file containing the mapping of col_name to loom file "path"')
parser.add_argument("--colname", default=None, type=str,
                    help='column name to use for extraction from samples file')
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

L.info("Reading sample.table")
sample_table = pd.read_csv(args.samples, sep="\t")

samples = dict(zip(sample_table[args.colname], sample_table.path))

# print(samples)


out_file = os.path.join(args.outdir, "matrix.loom")

# from https://linnarssonlab.org/loompy/cookbook/index.html

j = 1

with loompy.new(out_file) as dsout:  # Create a new, empty, loom file

    for seq_id in samples.keys():

        input_loom = samples[seq_id]

        L.info("Working on sample " + str(j) + ": " + input_loom)

        target_cells = cell_table["barcode"][cell_table[args.colname] == seq_id].values

        # add the sequencing identifier to the barcode

        target_cells = [x + "-" + seq_id for x in target_cells]
        # print(target_cells[0:10])

        L.info("Extracting " + str(len(target_cells)) + " cells")


        with loompy.connect(input_loom) as ds:

            # velocyto.py mangles the barcodes!
            # from AAAAAA-1
            # to   sample:AAAAAx


            # print(ds.ca['CellID'][0:10])

            # this is not good.. input files are being modified!
            if ":" in ds.ca['CellID'][0]:
                ds.ca['CellID'] = [x.split(":")[1].replace("x","-1-") + x.split(":")[0]
                                   for x in ds.ca['CellID']]

            mask = np.in1d(ds.ca['CellID'], target_cells)

            if np.any(mask) == False:
                raise ValueError("No matching cells found")

            for (ix, selection, view) in ds.scan(items=mask, axis=1, key="Accession"):
                dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

        j += 1

L.info("complete")
