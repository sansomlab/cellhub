
import numpy as np
import pandas as pd
import dandelion as ddl
import sys
import os
import logging
import argparse
from datetime import date
import glob
import anndata as ad
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("tcr_dandelion.py")


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--contig_path", default="", type=str,
                    help="path to the contig annotations")



parser.add_argument("--outfile_prefix",default=1, type=str,
                    help="prefix for output objects")

args = parser.parse_args()

L.info("args:")
print(args)

# ########################################################################### #
# ############################### Script    ################################# #
# ########################################################################### #


L.info("Running dandelion")

# Process each dandelion tsv file: 
if os.path.exists(args.contig_path):
    # Read in the output from dandelion preprocessing.
    # args.contig_path looks like this:
    # "dandelion/<sample_name>/dandelion/all_contig_dandelion.tsv"
    vdj = ddl.read_10x_airr(args.contig_path)

    # Add additional useful columns to metadata
    # duplicate_count represents number of umi
    ddl.update_metadata(
            vdj, retrieve="umi_count", retrieve_mode="split and merge"
    )

    # consensus_count represents number of reads
    ddl.update_metadata(
            vdj, retrieve="consensus_count", retrieve_mode="split and merge"
    )

    # Save an unfiltered version for more refined downstream processing
  
    # Save the metafile for downstream processing and attaching to the anndata object
    # vdj.metadata.to_csv( "dandeloin_meta.tsv", index=True, index_label="barcode", sep="\t")
    vdj.metadata.to_csv( args.outfile_prefix + "_meta_unfiltered.tsv", index=True, index_label="barcode", sep="\t")
    
    # Save the dandelion object for downstream processing and tcr visualisations
    vdj.write_h5ddl(args.outfile_prefix + "_unfiltered.h5ddl", complib="bzip2")

    # Perform filtering and clonotype calling for the filtered outputs
    vdj_filtered = ddl.pp.filter_contigs(vdj, 
        productive_only=True, 
        library_type="tr-ab",
        keep_highest_umi=True,
        umi_foldchange_cutoff=2,
        filter_extra_vdj_chains=True,
        filter_extra_vj_chains=False
        )


    ddl.tl.find_clones(vdj_filtered,identity=1, key="junction_aa")

    # Save the filtered metafile for downstream processing and attaching to the anndata object
    vdj_filtered.metadata.to_csv( args.outfile_prefix + "_meta_filtered.tsv", index=True, index_label="barcode", sep="\t")
    
    # Save the dandelion object for downstream processing and tcr visualisations
    vdj_filtered.write_h5ddl(args.outfile_prefix + "_filtered.h5ddl", complib="bzip2")

else: 
    raise ValueError("The contig path does not exist: " + args.contig_path)    


L.info("complete")