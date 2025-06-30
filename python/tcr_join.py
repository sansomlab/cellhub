
import numpy as np
import pandas as pd
import dandelion as ddl
from graph_tool.all import sfdp_layout
import sys
import os
import logging
import argparse
import scanpy as sc
from datetime import date
import glob
import anndata as ad
import matplotlib as mpl
# import psutil
mpl.rcParams.update(mpl.rcParamsDefault)


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("tcr_join.py")

# Attempt to retrieve SLURM_JOB_ID
slurm_job_id = os.getenv('SLURM_JOB_ID', 'Not set')
L.info(f'tcr_join.py SLURM_JOB_ID: {slurm_job_id}')
# print(f'tcr_join.py SLURM_JOB_ID: {slurm_job_id}')

# Attempt to get and print the 'PWD' environment variable
current_directory = os.environ.get('PWD')
L.info("The value of PWD after submission is: %s", current_directory)
# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--dandelion_path", default="", type=str,
                    help="path to the dandelion object")



parser.add_argument("--gex_path",default="", type=str,
                    help="path to the gex object")

                    
parser.add_argument("--cores_to_use",default=1, type=int,
                    help="Number of cores to use for parallel processing")

parser.add_argument("--outfile_prefix",default="", type=str,
                    help="The file prefix for the output files")

                    

args = parser.parse_args()

L.info("args: %s", args)



# ########################################################################### #
# ############################### Script    ################################# #
# ########################################################################### #


L.info("Running dandelion to join objects")



# # Get the number of logical CPUs in the system
# L.info("**********Logical CPUs:", psutil.cpu_count())

# # Get the total physical memory available
# L.info("**********Total memory:", psutil.virtual_memory().total)

# # Get the memory usage in percentage
# L.info("**********Memory usage (%):", psutil.virtual_memory().percent)

def subset_to_dandelion(adata, vdj):
    """
    Subsets an anndata object to only include those in a dandelion object.

    Parameters:
    adata: Anndata object.
    vdj:Dandelion object.

    Returns:
    a subsetted anndata object (adata_sub).
    """
    # Print initial number of indices
    print(f"Initial number of indices in adata: {adata.n_obs}")
    print(f"Initial number of indices in vdj: {len(vdj.data.cell_id)}")
    
    # Find the intersection of the indices
    common_indices = adata.obs_names.intersection(vdj.data.cell_id)
    print(f"Number of common indices: {len(common_indices)}")
    
    # Subset the anndata objects
    adata_sub = adata[common_indices].copy()
    print(f"New number of indices in adata: {adata_sub.n_obs}")
    
    return adata_sub

if os.path.exists(tmp_directory_name):
    L.info("Temporary directory exists")
else:
    L.info("Temporary directory does not exist")



# Process each dandelion object: 
if os.path.exists(args.dandelion_path):
    
    # Read in the dandelion object
    vdj = ddl.read_h5ddl(args.dandelion_path)

    # Filter dandelion object if specified


    if args.outfile_prefix == "dandelion_gex_merged_filtered":
        L.info("Filtering")
        vdj = ddl.pp.filter_contigs(vdj, 
            filter_contig=True, #Default
            library_type="tr-ab",
            filter_poorqualitycontig=False, #Default
            keep_highest_umi=True, #Default
            umi_foldchange_cutoff=2, #Default
            filter_extra_vdj_chains=True, #Default
            filter_extra_vj_chains=False, #Default
            productive_only=True #Default
        )

        vdj.data.shape 
        vdj.metadata.shape 
        ddl.update_metadata(vdj, retrieve="umi_count", retrieve_mode="split and merge" )
        vdj.data.shape
        vdj.metadata.shape
        vdj.write_h5ddl("vdj_prod_single.h5ddl", complib="bzip2")


    # Find clones
    ddl.tl.find_clones(vdj, identity=1, key="junction_aa", by_alleles=False, key_added="clone_id")

    vdj.write_h5ddl("vdj_prod_single_networked.h5ddl", complib="bzip2")



# COMMENT/UNCOMMMENT THIS SECTION TO ADD GEX DATA TO VDJ DATA - we generally don't do this anymore.
    if not os.path.exists(args.gex_path):
        raise ValueError("The GEX path does not exist: " + args.gex_path)
    L.info("Reading in GEX file")
    # Read in the GEX object
    adata = sc.read(args.gex_path)

    L.info("Subsetting GEX file to our TCR data")
    # subset adata to only cells with the same barcodes as vdj
    adata = subset_to_dandelion(adata, vdj)

    
    L.info("Saving subsetted GEX file for downstream analysis")
    # get the filename without the parent directories or extension of the args.gex_path variable


    # Extract the base filename with extension
    gex_filename = os.path.basename(args.gex_path)  # Results in "example.txt"

    # Split the extension and get the filename without extension
    gex_name = os.path.splitext(gex_filename)[0] 


    adata.write(args.outfile_prefix + "_" + gex_name + "_subset.h5ad")
    L.info("Finished saving subsetted GEX file for downstream analysis")

    L.info("Transferring vdj data to adata")
    ddl.tl.transfer(adata, vdj)
    L.info("Finished transferring vdj data to adata")


    L.info("Saving subsetted GEX file with TCR info for downstream analysis")
    #adata.write(args.outfile_prefix + "_subset_transfered.h5ad")
    adata.write(args.outfile_prefix + "_" + gex_name + "_subset_transfered.h5ad")
# END UNCOMMENT THIS SECTION TO ADD GEX DATA TO VDJ DATA


else: 
    raise ValueError("The dandelion path does not exist: " + args.dandelion_path)    


L.info("Complete")