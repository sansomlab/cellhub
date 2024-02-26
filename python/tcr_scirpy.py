
import numpy as np
import pandas as pd
#import scanpy as sc
import scirpy as ir
import sys
import os
import logging
import argparse
from datetime import date
import glob
import anndata as ad

#from cycler import cycler
#from matplotlib import cm as mpl_cm
#from matplotlib import pyplot as plt


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


parser.add_argument("--outfile_prefix",default=1, type=str,
                    help="name of the output metadata file")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################### Script.. ################################# #
# ########################################################################### #


L.info("Running scirpy")

# Step 1: Process each TCR csv file: 
if os.path.exists(args.contig_path):
    # Read the CSV file with ir.io.read_10x_vdj
    adata = ir.io.read_10x_vdj(args.contig_path)
else: 
    raise ValueError("The contig path does not exist: " + args.contig_path)    


# Step 2: Run chain quality control, adding to metadata
# Add chain primary and secondary VD/VDJ based on expression level,
# filter non-productive chains and chains without a valid CDR3,
# and add multichain label
# This is only available in a higher version of scirpy, which  requires python 3.9
##ir.pp.index_chains(adata)
# Add receptor_type, receptor_subtype and chain_pairing information
ir.tl.chain_qc(adata)


# Step 3: Add 'predicted_TCR_identity' column to adata 
def predict_tcr_identity(row):
    if row['IR_VJ_1_v_call'] == 'TRAV10' and row['IR_VJ_1_j_call'] == 'TRAJ18':
        return 'iNKT'
    elif row['IR_VJ_1_v_call'] == 'TRAV1-2' and row['IR_VJ_1_j_call'] in ['TRAJ12', 'TRAJ20', 'TRAJ33']:
        return 'MAIT'
    else:
        return 'NA'

adata.obs['predicted_TCR_identity'] = adata.obs.apply(predict_tcr_identity, axis=1) 

# Step 4: Add 'predicted_TCR_identity' column to adata
adata.obs['barcode'] = adata.obs.index

# ########################################################################### #
# ################################### save and exit ######################### #
# ########################################################################### #

L.info("saving scirpy metadata and AIRR data ")

# adata.obs.to_csv(args.outfile.replace(".sentinel", ".tsv.gz"), sep="\t", index=False)
adata.obs.to_csv(args.outfile_prefix + ".tsv.gz", sep="\t", index=False)
ir.io.write_airr(adata, args.outfile_prefix + "_AIRR.tsv")

L.info("complete")