'''
===================
Pipeline cellbender
===================

Overview
========
This pipeline uses CellBender to remove ambient UMI counts.

Configuration 
-------------
This pipeline is normally run in a seperate e.g. "cellhub_cellbender" directory so that the downstream results can be compared with those from cellranger.

The pipeline requires a configured :file:`pipeline_cellbender.yml` file. A default version can be obtained by executing: ::

   cellhub cellbender config


Input
-----

The location of the cellhub folder containing the cellranger results that will be used as the input for CellBender is be specified in the "pipeline_cellbender.yml" configuration file. Typically the user will have two parallel "cellhub" instances, e.g.:  

#. "cellhub" <- containing a first cellhub run based on the Cellranger counts (counts registered with "cellhub cellranger_multi make useCounts").

#. "cellhub_cellbender" <- containing a second cellhub run using CellBender to correct the Cellranger counts from the first run (counts registered with "cellhub cellbender make useCounts").


Running the pipeline
--------------------

It is recommended to run the cellbender task on a gpu queue.

On the University of Oxford's BMRC cluster, this can be achieved with e.g. ::

    cellhub cellbender make cellbender -v5 -p 200 --cluster-queue=short.qg --cluster-options "-l gpu=1,gputype=p100"


Pipeline output
===============

The pipeline registers cleaned CellBender h5 files on the local cellhub API. Currently this format is not fully compatible with the 10x h5 format. To work around this a custom loader is used, see the :doc:`cellhub.tasks.cellbender module documentation <../tasks/cellbender>` for more details.

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

import cellhub.tasks.parameters as chparam
import cellhub.tasks.TASK as TASK
import cellhub.tasks.api as api

import scanpy as sc

# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = chparam.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(chparam.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]

# ------------------------------ Pipeline Tasks ------------------------------ #


@follows(mkdir("cellbender.dir"))
@transform(glob.glob(os.path.join(PARAMS["cellhub_location"],
                                  "api/cellranger.multi/",
                                  "counts/unfiltered/*/h5/",
                                  "data.h5")),
           formatter(".*/data.h5"),
           r"cellbender.dir/{subdir[0][1]}/cellbender.sentinel")
def cellbender(infile, outfile):
    '''
    This task will run the CellBender comand.
    Please visit cellbender.readthedocs.io for further details.
    '''
    
    spec, SPEC = TASK.get_vars(infile, outfile, PARAMS)

    job_threads, job_memory, r_memory = TASK.get_resources(
        memory=PARAMS["resources_memory"], cpu=PARAMS["resources_cpu"],
        PARAMS=PARAMS)
        
    sample = str(os.path.basename(Path(infile).parents[1]))
    
    if PARAMS["cellbender_cuda"]:
        cuda_stat = "--cuda"
    else:
        cuda_stat = ""
        
    sample_key = "samples_" + sample #
   
    expected_cells = PARAMS[sample_key]["expected_cells"]
    total_droplets = PARAMS[sample_key]["total_droplets_included"]
    
    out_file = outfile.replace(".sentinel", ".h5")
    
    # Formulate and run statement
    statement = '''cellbender remove-background
                 --input=%(infile)s
                 --output=%(out_file)s 
                 --model=%(cellbender_model)s
                 %(cuda_stat)s
                 --expected-cells=%(expected_cells)s
                 --total-droplets-included=%(total_droplets)s
                 --fpr=%(cellbender_fpr)s
                 --epochs=%(cellbender_epochs)s
                 --learning-rate=%(cellbender_learning_rate)s
                 --low-count-threshold=%(cellbender_low_count_threshold)s
                 &> %(log_file)s
              ''' % dict(PARAMS, **SPEC, **locals())
    
    #print(statement)        
    P.run(statement)

    IOTools.touch_file(outfile)



@transform(cellbender,
           regex(r"cellbender.dir/(.*)/cellbender.sentinel"),
           r"cellbender.dir/\1/register.h5.sentinel")
def h5API(infile, outfile):
    '''
    Put the h5 files on the API

    Inputs:

        The input cellbender.dir folder layout is:

        unfiltered "outs": ::

            library_id/cellbender.h5

        filtered "outs": ::

            library_id/cellbender_filtered.h5

    '''
    x = api.api("cellbender")

    out_dir = os.path.dirname(outfile)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    library_id = infile.split("/")[1]

    h5_template = {"h5": {"path":"path/to/barcodes.tsv",
                          "format": "h5",
                          "link_name": "data.h5",
                          "description": "Cellbender h5 count file"}
                     }

    # 1. deal with unfiltered count data
    h5_location = os.path.join("cellbender.dir", library_id,
                                "cellbender.h5")

    h5_x = h5_template.copy()
    h5_x["h5"]["path"] = h5_location

    x.define_dataset(analysis_name="counts",
                     data_subset="unfiltered",
                     data_id=library_id,
                     data_format="h5",
                     file_set=h5_x,
                     analysis_description="Cellbender h5 file")


    x.register_dataset()

    # 2. deal with per sample libraries
    h5_location = os.path.join("cellbender.dir",
                               library_id,
                               "cellbender_filtered.h5")

    h5_x = h5_template.copy()
    h5_x["h5"]["path"] = h5_location

    x.define_dataset(analysis_name="counts",
                        data_subset="filtered",
                        data_id=library_id,
                        data_format="h5",
                        file_set=h5_x,
                        analysis_description="Cellbender h5 file")

    x.register_dataset()

    IOTools.touch_file(outfile)



@transform(cellbender,
           regex(r"cellbender.dir/(.*)/cellbender.sentinel"),
           r"cellbender.dir/\1/mtx.sentinel")
def mtx(infile, outfile):
    '''
        Convert cellbender h5 to mtx format
    '''

    spec, SPEC = TASK.get_vars(infile, outfile, PARAMS)

    job_threads, job_memory, r_memory = TASK.get_resources(
        memory=PARAMS["resources_memory"], cpu=PARAMS["resources_cpu"],
        PARAMS=PARAMS)
        
    statements = []
    
    to_process = {"filtered": "cellbender_filtered.h5",
                  "unfiltered": "cellbender.h5"}
    
    for type, path in to_process.items():
    
        h5 = os.path.join(spec.outdir, path)    
        mtx_dir = os.path.join(spec.outdir, type)
        
        log_name = spec.log_file.replace("mtx", "mtx." + type)
    
        # Formulate and run statement
        stat = '''python %(cellhub_code_dir)s/python/cellbender_export_mtx.py
                       --cellbender_h5=%(h5)s
                       --mtx_dir=%(mtx_dir)s
                     &> %(log_name)s
                    ''' % dict(PARAMS, **SPEC, **locals())
                    
        statements.append(stat)
    
    P.run(statements)

    IOTools.touch_file(outfile)


@transform(mtx,
           regex(r"cellbender.dir/(.*)/mtx.sentinel"),
           r"cellbender.dir/\1/register.mtx.sentinel")
def mtxAPI(infile, outfile):
    '''
    Put the mtx files on the API

    Inputs:

        The input cellbender.dir folder layout is:

        unfiltered "outs": ::

            library_id/cellbender.h5

        filtered "outs": ::

            library_id/cellbender_filtered.h5

    '''
    x = api.api("cellbender")

    out_dir = os.path.dirname(outfile)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    library_id = infile.split("/")[1]

    mtx_template = {"barcodes": {"path":"path/to/barcodes.tsv",
                                 "format": "tsv",
                                 "description": "cell barcode file"},
                    "features": {"path":"path/to/features.tsv",
                                  "format": "tsv",
                                  "description": "features file"},
                     "matrix": {"path":"path/to/matrix.mtx",
                                 "format": "market-matrix",
                                 "description": "Market matrix file"}
                     }
    
    to_register = {"unfiltered": os.path.join(out_dir, "unfiltered"),
                   "filtered": os.path.join(out_dir, "filtered")}
    
    for subset, mtx_loc in to_register.items():

        if os.path.exists(mtx_loc):

            mtx_x = mtx_template.copy()
            mtx_x["barcodes"]["path"] = os.path.join(mtx_loc, "barcodes.tsv.gz")
            mtx_x["features"]["path"] = os.path.join(mtx_loc, "features.tsv.gz")
            mtx_x["matrix"]["path"] =  os.path.join(mtx_loc, "matrix.mtx.gz")

            x.define_dataset(analysis_name="counts",
                             data_subset=subset,
                             data_id=library_id,
                             data_format="mtx",
                             file_set=mtx_x,
                             analysis_description="Cellbender count GEX output")

            x.register_dataset()

    IOTools.touch_file(outfile)
    

# -------------------------- Generic pipeline tasks -------------------------- #


@follows(h5API, mtxAPI)
def full():
    '''
    Run the full pipeline.
    '''
    pass
    
    
@follows(mtxAPI, h5API)
@files(None,"use.cellbender.sentinel")
def useCounts(infile, outfile):
    '''
        Set the cellbender counts as the source for downstream analysis.
        This task is not run by default.
    '''
    
    if os.path.exists("api/counts"):
        raise ValueError("Counts have already been registered to the API")

    else:
        os.symlink("cellbender/counts", "api/counts")


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
