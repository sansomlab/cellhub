
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
L = logging.getLogger("dandelion_join.py")

# Attempt to retrieve SLURM_JOB_ID
slurm_job_id = os.getenv('SLURM_JOB_ID', 'Not set')
L.info(f'dandelion_join.py SLURM_JOB_ID: {slurm_job_id}')
# print(f'dandelion_join.py SLURM_JOB_ID: {slurm_job_id}')


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--dandelion_input_files", default="", type=str,
                    help="paths to the dandelion objects to be concatenated")



parser.add_argument("--output_file",default="", type=str,
                    help="The output sentinel filename")

                    

args = parser.parse_args()

L.info("args: %s", args)



# ########################################################################### #
# ############################### Script    ################################# #
# ########################################################################### #


L.info("Running dandelion to concatenate dandelion objects")

# Convert the dandelion_input_files string, where each item is seperated by a semicolon, to a list
input_files_list = args.dandelion_input_files.split(",")
output_file = args.output_file

# # Get the number of logical CPUs in the system
# L.info("**********Logical CPUs:", psutil.cpu_count())

# # Get the total physical memory available
# L.info("**********Total memory:", psutil.virtual_memory().total)

# # Get the memory usage in percentage
# L.info("**********Memory usage (%):", psutil.virtual_memory().percent)

if(len(input_files_list) == 0):
    raise ValueError("No dandelion files found")
else:
    vdj_list = []
    for file_path in input_files_list:
        print(file_path)
        vdj = ddl.read_h5ddl(file_path)
        vdj_list.append(vdj)
    
    if(len(vdj_list) == 1):
        vdj = vdj_list[0]
    else:
        vdj = ddl.concat(vdj_list)


    #########
    # We have to re-calculate the additional metadata as these do get concatenated, only the data gets concatenated

    # Add additional useful columns to metadata
    # duplicate_count / umi_count represents number of umi 
    ddl.update_metadata(
            vdj, retrieve="umi_count", retrieve_mode="split and merge"
    )

    # consensus_count represents number of reads
    ddl.update_metadata(
            vdj, retrieve="consensus_count", retrieve_mode="split and merge"
    )


    # Save the metafile and ddl object for downstream processing and attaching to the anndata object
    vdj.metadata.to_csv( output_file.replace(".sentinel", "_unfiltered_meta.tsv"), index=True, index_label="barcode", sep="\t")
    vdj.write_h5ddl(output_file.replace(".sentinel", "_unfiltered.h5ddl"), complib="bzip2")

    # Also save the AIRR file for downstream processing (needs testing, moved from TCR match)
    vdj.write_airr(output_file.replace(".sentinel", "_unfiltered_AIRR.tsv"))    


L.info("Complete")