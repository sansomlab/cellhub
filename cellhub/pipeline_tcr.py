'''
======================
pipeline_tcr.py
======================


Overview
========

This pipeline performs the following functions:

* Reproduce combat?
* Scirpy?

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` for general
information on how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured:file:`pipeline_cellranger.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_cellranger.py config


Inputs
------

In addition to the "pipeline_tcr.yml" file, the pipeline requires two inputs: 

# api/cellranger/ .... tcr (vdj-t) ... 

Dependencies
------------

This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* cellranger: https://support.10xgenomics.com/single-cell-gene-expression/


Pipeline logic
--------------

The pipeline is designed to:

* 

Pipeline output
---------------

The pipeline returns:

* tsv tables that map the cell barcode to tcr-related statistics and annotations.

* e.g. santised Cellranger vdj output/ scirpy output/IMGT output....

Code
====

'''

from ruffus import *
from pathlib import Path
import sys
import os
import glob
import sqlite3
import yaml
import  csv

import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools

import pandas as pd
import numpy as np

# import local pipeline utility functions
import cellhub.tasks as T
import cellhub.tasks.cellranger as cellranger
import cellhub.tasks.samples as samples

# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = T.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(T.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]


# ----------------------- Read in the samples set -------------------------- #

runCount, runTCR, runBCR = False, False, False

# Only do this when tasks are being executed.
if len(sys.argv) > 1:
    if sys.argv[1] == "make":
        
        S = samples.samples(sample_tsv = PARAMS["sample_table"],
                    library_tsv = PARAMS["library_table"])
        
        if any([x in S.known_feature_library_types  
                for x in S.feature_types]): 
            runCount = True

        if "VDJ-T" in S.feature_types: runTCR = True 
        if "VDJ-B" in S.feature_types: runBCR = True 

# ---------------------------- Pipeline tasks ------------------------------- #

# ########################################################################### #
# ###########################  TCR Analysis  ############################## #
# ########################################################################### #


# 1. recreation of basic stats ala COMBAT with e.g. scirpy or manual script (TCR QC)


@follows(mkdir("cell.qc.dir"))
@transform(glob.glob("api/tcr/path/libraryX/*.cellranger.per.library.tables"),  ## FIX PATH
           regex(r".*/.*/.*/(.*)/mtx/matrix.mtx.gz"), ## 
           r"scirpy.dir/\1/scirpy.sentinel")
def scirpy(infile, outfile):

    t = T.setup(infile, outfile, PARAMS,
                memory=PARAMS["cellranger_localmem"],
                cpu=PARAMS["cellranger_localcores"])

    #t.outdir == outdir = os.path.dirname(outfile)
    #t.resources == {"job_threads": cores; "job_memory": G}
  

    statement = '''python %(cellhub_code_dir)s/python/tcr_scirpy.py
                   -- inputs?
                   -- params?
                   -- outputs %(outdir)s
                    &> %(log_file)s
                 ''' % dict(PARAMS, 
                            **t.var, 
                            **locals())

    P.run(statement, **t.resources)
    IOTools.touch_file(outfile)


# # In this section the pipeline processes the gene expression (GEX) and antibody
# # capture, i.e. antibody derived tag (ADT) information.

# def count_jobs():

#     if not os.path.exists("cellranger.count.dir"):
#         os.mkdir("cellranger.count.dir")
    
#     for lib in S.feature_barcode_libraries():
    
#         csv_path = os.path.join("cellranger.count.dir", lib + ".csv")
    
#         if not os.path.exists(csv_path):
#             S.write_csv(lib, csv_path)
    
#         yield(csv_path, os.path.join("cellranger.count.dir",
#                                  lib + ".sentinel"))
    
# @active_if(runCount)      
# @files(count_jobs)
# def count(infile, outfile):
#     '''
#     Execute the cellranger count pipeline
#     '''
    
#     t = T.setup(infile, outfile, PARAMS,
#                 memory=PARAMS["cellranger_localmem"],
#                 cpu=PARAMS["cellranger_localcores"])

#     this_library_id = os.path.basename(infile)[:-len(".csv")]

#     library_parameters = S.samples[this_library_id]

#     # provide references for the present feature types
#     lib_types = S.lib_types(this_library_id)
#     transcriptome, feature_ref  = "", ""
    
#     if "Gene Expression" in lib_types:
#         transcriptome = "--transcriptome=" + PARAMS["gex_reference"]
    
#     if "Antibody Capture" in lib_types:
#         feature_ref =  "--feature-ref" + PARAMS["feature_reference"]

#     # add read trimming if specified
#     r1len, r2len = "", ""

#     if PARAMS["gex_r1-length"] != "false":
#         r1len = PARAMS["gex_r1-length"]
    
#     if PARAMS["gex_r2-length"] != "false":
#         r1len = PARAMS["gex_r2-length"]
 
#     # deal with flags
#     nosecondary, nobam, includeintrons = "", "", ""
    
#     if PARAMS["cellranger_nosecondary"]:
#         nosecondary = "--nosecondary"
#     if PARAMS["cellranger_no-bam"]:
#         nobam = "--no-bam"
#     if PARAMS["gex_include-introns"]:
#         includeintrons = "--include-introns=true"
 
#     statement = '''cd cellranger.count.dir;
#                     cellranger count
# 	    	        --id %(this_library_id)s
#                     %(transcriptome)s
#                     %(feature_ref)s
#                     --libraries=../%(infile)s
# 		            --nopreflight
#                     --disable-ui
#                     --expect-cells=%(expect_cells)s
#                     --chemistry=%(chemistry)s
#                     %(nosecondary)s
#                     %(nobam)s
#                     --localcores=%(cellranger_localcores)s
#                     --localmem=%(cellranger_localmem)s
#                     %(includeintrons)s
#                     %(r1len)s %(r2len)s
#                     &> ../%(log_file)s
#                  ''' % dict(PARAMS, 
#                             **library_parameters,
#                             **t.var, 
#                             **locals())

#     P.run(statement, **t.resources)
#     IOTools.touch_file(outfile)




# ---------------------------< Pipeline targets >------------------------------ #

@follows(count, 
         mtxAPI, h5API, 
         registerTCR, registerBCR,
         registerMergedBCR, registerMergedTCR)
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
