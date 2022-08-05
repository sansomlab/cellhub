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



# <---------------------------- Parse arguments ----------------------------> #

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object"),
parser.add_argument("--conserved", action='store_true',
                    help="find conserved markers")   
parser.add_argument("--conserved_factor", default="None", type=str,
                    help="name of the conserved factor"), 
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output metadata file ")
args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# <--------------------------------------------------------------------------> #

adata = ad.read_h5ad(args.source_anndata, backed='r')

L.info("Saving the metadata")
adata.obs.to_csv(args.outfile, sep="\t", index=False)    

outdir = os.path.dirname(args.outfile)

if args.conserved:

    L.info('Writing out the levels of the conserved factor')

    levels = [x for x in adata.obs[args.conserved_factor].cat.categories]
    print(levels)

    levels_file = os.path.join(outdir,
                               args.conserved_factor + ".levels")

    with open(levels_file,"w") as cons_levels:
        for level in levels:
           cons_levels.write(level + "\n")

# <------------------------------- functions --------------------------------> #
                   

L.info("Complete")
