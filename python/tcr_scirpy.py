
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

#parser.add_argument("--library_id", default="umap.tsv.gz", type=str,
#                    help="a comma seperated list of umaps to add")

parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################### Script.. ################################# #
# ########################################################################### #


L.info("Running scirpy")



# Step 1: Get contigs from the cellhub directory

# df_libraries = pd.read_csv("libraries.tsv", sep='\t')
# filtered_libraries = df_libraries[df_libraries['feature_type'] == 'VDJ-T']

# Step 1: Process each filtered library: 

    
if os.path.exists(args.contig_path):
        # Read CSV file into a pandas DataFrame
        # df = pd.read_csv(file_path, compression='gzip')
        
        # # Check if 'barcode' column exists, otherwise rename 'barcode_id': 
        # # need to check if cellranger multi generates a df with column barcode or barcode_id
        # if 'barcode' not in df.columns and 'barcode_id' in df.columns:
        #     df.rename(columns={'barcode_id': 'barcode'}, inplace=True)
            
        
        # # Save the modified DataFrame as CSV
        # filtered_csv_path = os.path.expanduser(f"filtered_{library_id}.csv")
        # df.to_csv(filtered_csv_path, index=False)
        
        # Read the CSV file with ir.io.read_10x_vdj
    adata = ir.io.read_10x_vdj(args.contig_path)
else: 
    raise ValueError("The contig path does not exist: " + args.contig_path)    


# Step 2: Run chain quality control
ir.tl.chain_qc(adata)

# Step 5: Add 'predicted_TCR_identity' column to adata 
def predict_tcr_identity(row):
    if row['IR_VJ_1_v_call'] == 'TRAV10' and row['IR_VJ_1_j_call'] == 'TRAJ18':
        return 'iNKT'
    elif row['IR_VJ_1_v_call'] == 'TRAV1-2' and row['IR_VJ_1_j_call'] in ['TRAJ12', 'TRAJ20', 'TRAJ33']:
        return 'MAIT'
    else:
        return 'NA'

adata.obs['predicted_TCR_identity'] = adata.obs.apply(predict_tcr_identity, axis=1) 

adata.obs['barcode'] = adata.obs.index

# ########################################################################### #
# ################################### save and exit ######################### #
# ########################################################################### #

L.info("saving the object ")

# with open() as ...:
#   for cell in exp:
#       write(result)

adata.obs.to_csv(args.outfile, sep="\t", index=False)

L.info("all done")