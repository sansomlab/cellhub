'''
================
Pipeline ADT dsb normalization
================

Overview
========
This pipeline performs the following steps:

* Calculates per-cell total_UMI for the unfiltered GEX and ADT modalities
* Pick up UMI depth to distinguish empty-drop from cell-containing barcodes
* Perform dsb normalization


Configuration
-------------
The pipeline requires a configured :file:`pipeline_adt_norm.yml` file. Default configuration files can be generated by executing: ::

   python <srcdir>/pipeline_adt_norm.py config


Input files
-----------
./api/cellranger.multi/ADT/unfiltered/*/mtx/*.gz 
./api/cellranger.multi/GEX/unfiltered/*/mtx/*.gz 


Dependencies
------------
This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* R dependencies required in the r scripts


Pipeline output
===============
The pipeline returns:
* adt_dbs.dir folder with per-input umidepth.tsv.gz table and a folder per sample with the normalized market matrices

Code
====

'''


from ruffus import *
from ruffus.combinatorics import *
import sys
import os
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools
from pathlib import Path
import pandas as pd
import glob

import cellhub.tasks.control as C
import cellhub.tasks.api as api

# Override function to collect config files
P.control.write_config_files = C.write_config_files


# -------------------------- < parse parameters > --------------------------- #

# load options from the yml file
parameter_file = C.get_parameter_file(__file__, __name__)
PARAMS = P.get_parameters(parameter_file)

# Set the location of the cellhub code directory
if "code_dir" not in PARAMS.keys():
    PARAMS["code_dir"] = Path(__file__).parents[1]
else:
    if PARAMS["code_dir"] != Path(__file__).parents[1]:
        raise ValueError("Could not set the location of "
                         "the pipeline code directory")

# ----------------------- < pipeline configuration > ------------------------ #

# handle pipeline configuration
if len(sys.argv) > 1:
        if(sys.argv[1] == "config") and __name__ == "__main__":
                    sys.exit(P.main(sys.argv))


# ########################################################### #
# ######## Calculate seq depth distributions  ############### #
# ########################################################### #

@follows(mkdir("adt_dsb.dir"))
@transform(glob.glob("api/cellranger.multi/GEX/unfiltered/*/mtx/matrix.mtx.gz"),
           regex(r".*/.*/.*/.*/(.*)/mtx/matrix.mtx.gz"),
           r"adt_dsb.dir/\1/\1_gex.sentinel")
def gexdepth(infile, outfile):
    '''This task will run R/calculate_depth_dist.R,
    It will describe the UMI distribution and a split it among background and cell-containing barcodes
    - Calculate the number of UMI per barcode
    - Label barcodes as True/False based on whether they are part or not of the background
    '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Get cellranger directory and id
    library_name = os.path.basename(outfile)[:-len(".sentinel")]
    cellranger_dir = os.path.dirname(infile)

    # Other settings
    job_threads = PARAMS["resources_threads"]
    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    log_file = outfile.replace(".tsv.gz", ".log")

    out_file = outfile.replace(".sentinel", ".tsv.gz")

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/calculate_depth_dist.R
                 --cellranger_dir=%(cellranger_dir)s
                 --library_id=%(library_name)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(out_file)s
              '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)


@merge(gexdepth,
       "adt_dsb.dir/api_gex.sentinel")
def gexdepthAPI(infiles, outfile):
    '''
    Add the umi depth metrics results to the API
    '''

    file_set={}

    for lib in infiles:

        tsv_path = lib.replace(".sentinel",".tsv.gz")
        library_id = os.path.basename(tsv_path).replace("_gex", "")

        file_set[library_id] = {"path": tsv_path,
                                "description":"all barcodes gex umi depth table for library " +\
                                library_id,
                                "format":"tsv"}

    x = api.api("adt_norm")

    x.define_dataset(analysis_name="depth_metrics",
              data_subset="gex",
              file_set=file_set,
              analysis_description="per library tables of cell GEX depth")

    x.register_dataset()

    # Create sentinel file
    IOTools.touch_file(outfile)


@follows(mkdir("adt_dsb.dir"))
@transform(glob.glob("api/cellranger.multi/ADT/unfiltered/*/mtx/matrix.mtx.gz"),
           regex(r".*/.*/.*/.*/(.*)/mtx/matrix.mtx.gz"),
           r"adt_dsb.dir/\1/\1_adt.sentinel")
def adtdepth(infile, outfile):
    '''This task will run R/calculate_depth_dist.R,
    It will describe the UMI distribution and a split it among background and cell-containing barcodes
    - Calculate the number of UMI per barcode
    - Label barcodes as "background" or "cell" based on whether they are part or not of the background
    '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Get cellranger directory and id
    library_name = os.path.basename(outfile)[:-len(".sentinel")]
    cellranger_dir = os.path.dirname(infile)

    # Other settings
    job_threads = PARAMS["resources_threads"]
    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    log_file = outfile.replace(".tsv.gz", ".log")

    out_file = outfile.replace(".sentinel", ".tsv.gz")

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/calculate_depth_dist.R
                 --cellranger_dir=%(cellranger_dir)s
                 --library_id=%(library_name)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(out_file)s
              '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)


@merge(adtdepth,
       "adt_dsb.dir/api_adt.sentinel")
def adtdepthAPI(infiles, outfile):
    '''
    Add the umi depth metrics results to the API
    '''

    file_set={}

    for lib in infiles:

        tsv_path = lib.replace(".sentinel",".tsv.gz")
        library_id = os.path.basename(tsv_path).replace("_adt", "")

        file_set[library_id] = {"path": tsv_path,
                                "description":"all barcodes adt umi depth table for library " +\
                                library_id,
                                "format":"tsv"}

    x = api.api("adt_norm")

    x.define_dataset(analysis_name="depth_metrics",
              data_subset="adt",
              file_set=file_set,
              analysis_description="per library tables of cell ADT depth")

    x.register_dataset()

    # Create sentinel file
    IOTools.touch_file(outfile)



@follows(adtdepthAPI)
@transform(glob.glob("api/cellranger.multi/ADT/unfiltered/*/mtx/matrix.mtx.gz"),
           regex(r".*/.*/.*/.*/(.*)/mtx/matrix.mtx.gz"),
           r"adt_dsb.dir/\1/mtx/\1.sentinel")
def dsb_norm(infile, outfile):
    '''This task will run R/normalize_adt.R,
    It will read the infiltered ADT count matrix. If not user definition of background and cell containig 
    barcodes, then the automatic guess from the gex and adt get depth tasks will be used.
    - Calculate dsb normalized ADT expression matrix
    - Write market matrices per sample
    '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Get cellranger directory and id
    library_name = os.path.basename(outfile)[:-len(".sentinel")]
    cellranger_dir = os.path.dirname(infile)

    # Other settings
    job_threads = PARAMS["resources_threads"]
    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    gex_depth = "api/adt.norm/depth_metrics/gex/" + library_name + "_gex.tsv.gz"
    adt_depth = "api/adt.norm/depth_metrics/adt/" + library_name + "_adt.tsv.gz"
 
    #if (PARAMS["dsb_background_counts_min"]) : 
    #    bcmin = PARAMS["dsb_background_counts_min"]
    #    bcmax = PARAMS["dsb_background_counts_max"]
    #    bfmin = PARAMS["dsb_background_feats_min"]
    #    bfmax = PARAMS["dsb_background_feats_max"]
    #    ccmin = PARAMS["dsb_cell_counts_min"]
    #    ccmax = PARAMS["dsb_cell_counts_max"]
    #    cfmin = PARAMS["dsb_cell_feats_min"]
    #    cfmax = PARAMS["dsb_cell_feats_max"]

    log_file = outfile.replace(".tsv.gz", ".log")

    out_file = outfile.replace(".sentinel", ".mtx.gz")

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/normalize_adt.R
                 --cellranger_dir=%(cellranger_dir)s
                 --library_id=%(library_name)s
                 --gex_depth=%(gex_depth)s
                 --adt_depth=%(adt_depth)s
                 --bcmin=%(bcmin)s
                 --bcmax=%(bcmax)s
                 --bfmin=%(bfmin)s
                 --bfmax=%(bfmax)s
                 --ccmin=%(ccmin)s
                 --ccmax=%(ccmax)s
                 --cfmin=%(cfmin)s
                 --cfmax=%(cfmax)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(out_file)s
              '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)


@merge(adtdepth,
       "adt_dsb.dir/api_adt.sentinel")
def dsbAPI(infiles, outfile):
    '''
    Add the umi depth metrics results to the API
    '''

    file_set={}

    for libqc in infiles:

        tsv_path = libqc.replace(".sentinel",".tsv.gz")
        library_id = os.path.basename(tsv_path)

        file_set[library_id] = {"path": tsv_path,
                                "description":"all barcodes adt umi depth table for library " +\
                                library_id,
                                "format":"tsv"}

    x = api.api("adt_norm")

    x.define_dataset(analysis_name="qcmetrics",
              data_subset="unfiltered",
              file_set=file_set,
              analysis_description="per library tables of cell ADT depth")

    x.register_dataset()


# ---------------------------------------------------
# Generic pipeline tasks

@follows(mkdir("adt_dsb.dir"))
@files(None, "adt_dsb.dir/plot.sentinel")
def plot(infile, outfile):
    '''Draw the pipeline flowchart'''

    pipeline_printout_graph ( "adt_dsb.dir/pipeline_flowchart.svg",
                          "svg",
                          [full],
                          no_key_legend=True)

    pipeline_printout_graph ( "adt_dsb.dir/pipeline_flowchart.png",
                          "png",
                          [full],
                          no_key_legend=True)

    IOTools.touch_file(outfile)


@follows(gexdepthAPI, adtdepthAPI, dsb_norm, plot)
def full():
    '''
    Run the full pipeline.
    '''
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
