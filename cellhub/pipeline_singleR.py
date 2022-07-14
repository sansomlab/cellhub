"""
================
Pipeline singleR
================

Overview
========

This pipeline runs `singleR<https://bioconductor.org/packages/release/bioc/html/SingleR.html>`_ for cell prediction. Single R:

(1) runs at cell level (cells are scored independently)
(2) Uses a non-paramentric correlation test (i.e. monotonic transformations of 
    the test data have no effect).

Given these facts, in cellhub we run singleR on the raw counts upstream to 
(a) help with cell QC and (b) save time in the interpretation phase.

This pipeline operates on the ensembl_ids.

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline should be run in the cellhub directory.

To obtain a configuration file run "cellhub singleR config".


Inputs
------

1. Per-sample h5 files (from the cellhub API).

2. References for singleR  obtained via the R bioconductor 'celldex' library. 
   As downloading of the references is very slow, they need to 
   be manually downloaded and "stashed" as rds files in an appropriate location using the 
   R/scripts/singleR_stash_references.R scripts. This location is then specified in the
   yaml file.


Pipeline output
===============

The pipeline saves the singleR scores and predictions for each of the
specified references on the cellhub API.

Code
====

"""

import os
import sys
import gzip
from shutil import copyfile
import glob

from pathlib import Path
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools

import cellhub.tasks.control as C
import cellhub.tasks.TASK as TASK
import cellhub.tasks.fetch_cells as fetch_cells

import cellhub.tasks.api as api

# Override function to collect config files
P.control.write_config_files = C.write_config_files


# -------------------------- < parse parameters > --------------------------- #

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline_singleR.yml" % os.path.splitext(__file__)[0],
     "../pipeline_singleR.yml",
     "pipeline_singleR.yml"])

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]


# ----------------------- < pipeline configuration > ------------------------ #

# handle pipeline configuration
if len(sys.argv) > 1:
        if(sys.argv[1] == "config") and __name__ == "__main__":
                    sys.exit(P.main(sys.argv))

# ------------------------------ < functions > ------------------------------ #

# ########################################################################### #
# ############################# pipeline tasks ############################## #
# ########################################################################### #

def genSingleRjobs():
    '''
       generate the singleR jobs
    '''

    sample_h5s = glob.glob("api/cellranger.multi/counts/filtered/*/h5/sample_filtered_feature_bc_matrix.h5")

    references = [x.strip() for x in PARAMS["reference_data"].split(",")]

    for h5_path in sample_h5s:

        sample = os.path.basename(Path(h5_path).parents[1])

        for reference in references:

            yield [h5_path,
                   os.path.join("singleR.dir",
                                reference,
                                sample + ".sentinel")]


@files(genSingleRjobs)
def singleR(infile, outfile):
    '''
       Perform cell identity prediction with singleR.
    '''
    
    spec, SPEC = TASK.get_vars(infile, outfile, PARAMS)

    job_threads, job_memory, r_memory = TASK.get_resources(
        cpu=PARAMS["resources_cores"],
        memory=PARAMS["resources_memory"])

    reference = os.path.basename(
    Path(outfile).parents[0]).replace(".ref.dir","")

    sample = os.path.basename(outfile).split(".")[0]

    statement = '''Rscript %(cellhub_code_dir)s/R/scripts/singleR_run.R
                        --h5file=%(infile)s
                        --sample=%(sample)s
                        --reference=%(reference)s
                        --refstashdir=%(reference_stash_dir)s
                        --workers=%(resources_cores)s
                        --outdir=%(outdir)s
                        &> %(log_file)s
                        ''' % dict(PARAMS, **SPEC, **locals())

    P.run(statement)
    IOTools.touch_file(outfile)


@follows(singleR)
@files(None,
       "singleR.dir/out.dir/labels.sentinel")
def concatenate(infile, outfile):
    '''
       Concatenate the label predictions across all the samples.
    '''
    
    spec, SPEC = TASK.get_vars(infile, outfile, PARAMS)

    job_threads, job_memory, r_memory = TASK.get_resources(
        memory=PARAMS["resources_memory"])
    
    
    references = [x.strip() for x in PARAMS["reference_data"].split(",")]
    
    stats = []
    
    for reference in references:
        ref_path = os.path.join("singleR.dir", reference)

        outdir = os.path.join(spec.outdir, reference)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        label_glob = "*.labels.tsv.gz" 
        label_table = os.path.join(outdir, "labels.tsv.gz")
        
        statement = '''zcat %(ref_path)s/%(label_glob)s
                       | awk 'NR!=1 { while (/^barcode_id/) getline; } 1 {print}'  
                       | gzip -c > %(label_table)s
                    '''  % locals()
                    
        stats.append(statement)
        
        score_glob = "*.scores.tsv.gz" 
        score_table = os.path.join(outdir, "scores.tsv.gz")
        
        statement = '''zcat %(ref_path)s/%(score_glob)s
                       | awk 'NR!=1 { while (/^barcode_id/) getline; } 1 {print}'  
                       | gzip -c > %(score_table)s
                    '''  % locals()
                    
        stats.append(statement)
    
    P.run(stats)
    IOTools.touch_file(outfile)


@files(concatenate,
       "singleR.dir/api.sentinel")
def singleRapi(infiles, outfile):
    '''
        Add the singleR results to the cellhub API.
    '''
    references = [x.strip() for x in PARAMS["reference_data"].split(",")]
    
    for reference in references:
    
        file_set={}
           
        file_set["labels"] = {
            "path": os.path.join("singleR.dir", "out.dir",
            reference, "labels.tsv.gz"),
            "description":"single R label predictions",
            "format":"tsv"}
        
        file_set["scores"] = {
            "path": os.path.join("singleR.dir", "out.dir",
            reference, "scores.tsv.gz"),
            "description":"single R prediction scores",
            "format":"tsv"}
        
        x = api.api("singleR")

        x.define_dataset(analysis_name=reference,
              #data_subset=reference,
              file_set=file_set,
              analysis_description="Aggregated SingleR prediction results")

        x.register_dataset()


# ########################################################################### #
# ##################### full target: to run all tasks ####################### #
# ########################################################################### #


@follows(singleR)
def full():
    pass


# ------------------- < ***** end of pipeline **** > ------------------------ #


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
