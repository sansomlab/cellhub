#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import os
import harmonypy as hm
import pandas as pd
import matplotlib.pyplot as plt
import logging
import yaml

# ########################################################################### #
# ############### Set up the log and figure folder ########################## #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("run_lisi")

# ########################################################################### #
# ############################ Script arguments ############################# #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--task-yml", default="", type=str,
                    help="yml for this task")
args = parser.parse_args()


# Read YAML file
with open(args.task_yml, 'r') as stream:
    opt = yaml.safe_load(stream)

L.info("Running with options ---> %s", opt)

# ########################################################################### #
# ######################## Read input data ################################## #
# ########################################################################### #

comp_file = opt["comp_file"]

components = pd.read_csv(comp_file, sep = "\t")
components.set_index("barcode", inplace=True)

metadata = pd.read_csv(os.path.join(opt["matrixdir"], "metadata.tsv.gz"), sep = "\t")

L.info("Read input components and metadata")

# ########################################################################### #
# ################################# Run LISI ################################ #
# ########################################################################### #

select = ['barcode'] + [(opt["split_var"])]
metadata_lisi = metadata.loc[:, select]

L.info("Compute lisi")

lisi = hm.compute_lisi(components, metadata_lisi, [(opt["split_var"])])

lisi = pd.DataFrame(lisi)
lisi.columns = ["iLISI"]
lisi['barcode'] = metadata_lisi.loc[:,'barcode']

lisi.to_csv(os.path.join(opt["outdir"], "lisi.tsv.gz"),
            sep="\t", index=False, compression="gzip")

L.info("Done running and saving LISI")

figdir = os.path.join(opt["outdir"], "figures.dir")
if not os.path.exists(figdir):
    os.makedirs(figdir)

# make a simple histogram
lisi.plot.hist(grid=True, bins=20, rwidth=0.9,
                   color='#607c8e')
plt.xlabel('iLISI')
plt.ylabel('n cells')
plt.grid(axis='y', alpha=0.75)
plt.savefig(os.path.join(figdir, "histogram_lisi.pdf"))

L.info("Done plotting a simple histogram of LISI values")

L.info("Completed")
