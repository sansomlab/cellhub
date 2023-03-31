'''
==================
pipeline_dehash.py
==================


Overview
========

This pipeline dehashes cells which have multiplexed using hash tag oligos (HTOs).

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` for general information on how to use cgat
pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

Default configuration files can be generated by executing:
   cellhub celldb config

Code
====

'''

from ruffus import *

import sys
import os
import re
import sqlite3
import pandas as pd
import numpy as np
import glob
from pathlib import Path

import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools
import cgatcore.database as database

import cellhub.tasks as T
import cellhub.tasks.dehash as dehash

# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = T.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(T.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]

# ------------------------------ Pipeline Tasks ------------------------------ #

@transform("api/counts/filtered/*/mtx/matrix.mtx.gz",
           regex("api/counts/filtered/(.*)/mtx/matrix.mtx.gz"),
           r"dehash.dir/gmm.demux.dir/\1.gmm.demux.sentinel")
def gmmDemux(infile, outfile):
    '''
    Run gmmDemux
    '''
    
    t = T.setup(infile, outfile, PARAMS,
                memory = PARAMS["hto_memory"],
                make_outdir=False)

    library_id = os.path.basename(outfile)[:-len(".gmm.demux.sentinel")]

    input_mtx = os.path.dirname(infile)

    gmm_working_dir = os.path.join(os.path.dirname(outfile),
                                   library_id)

    if not os.path.exists(gmm_working_dir):
        os.makedirs(gmm_working_dir)

    if PARAMS["hto_per_library"] == True:
        #HTOs = "_".join([PARAMS["hto"], library_id])
        HTOs = PARAMS["hto_"+library_id]
    else:
        HTOs = PARAMS["hto_names"]

    if PARAMS["gmm_demux_per_library"] == True:
        
        threshold = PARAMS["gmm_demux_"+library_id]
    else:
        threshold = PARAMS['gmm_demux_threshold']

    statement = '''GMM-demux %(input_mtx)s
                             %(HTOs)s
                             --threshold %(threshold)s
                             --full %(gmm_working_dir)s/full
                             --simplified %(gmm_working_dir)s/simple
                             --output %(gmm_working_dir)s/SSD
                             --summary %(gmm_demux_ncells)s
                             --report %(gmm_working_dir)s/report
                &> %(log_file)s
                ''' % dict(PARAMS, **t.var, **locals())

    P.run(statement, **t.resources)

    # parse the output to something more sensible
    # output has columns
    # gmm_cluster_id, gmm_confidence, gmm_call, barcode_id, gmm_singlet
    hto_names = [x.strip() for x in HTOs.split(",")]

    results_dir = "dehash.dir/gmm.demux.dir/results.dir/"

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    dehash.parse_gmmdemux(gmm_working_dir,
                          results_dir,
                          library_id,
                          hto_names)

    IOTools.touch_file(outfile)


@merge(gmmDemux,
       "dehash.dir/gmm.demux.api.sentinel")
def gmmAPI(infiles, outfile):
    '''
    Register the GMM-Demux results on the API
    '''

    file_set={}

    for x in infiles:


        library_id = os.path.basename(x)[:-len(".gmm.demux.sentinel")]

        tsv_path = os.path.join(os.path.dirname(x),
                                "results.dir",
                                library_id + ".tsv.gz")

        file_set[library_id] = {"path": tsv_path,
                                "description":"Parsed GMM-demux results for " +\
                                library_id,
                                "format":"tsv"}

    x = T.api("dehash")

    x.define_dataset(analysis_name="gmm.demux",
              data_subset="filtered",
              file_set=file_set,
              analysis_description="per library tables of GMM-demux results")

    x.register_dataset()


@transform("api/counts/filtered/*/mtx/matrix.mtx.gz",
           regex("api/counts/filtered/(.*)/mtx/matrix.mtx.gz"),
           r"dehash.dir/demuxEM.dir/\1.demuxEM.hash.count.csv.sentinel")
def hashCountCSV(infile, outfile):
    '''
    Make the hash count csv table for demuxEM
    '''

    demuxEM_working_dir = os.path.dirname(outfile)
    if not os.path.exists(demuxEM_working_dir):
        os.makedirs(demuxEM_working_dir)

    import anndata as anndata
    import scanpy as sc
    import pandas as pd

    if PARAMS["hto_per_library"] == True:
        HTOs = "_".join([PARAMS["hto"], library_id])

    else:
        HTOs = PARAMS["hto_names"]

    ## move this to a helper script.
    ad = sc.read_10x_mtx(os.path.dirname(infile), gex_only=False)
    ad = ad.T
    ad = ad[HTOs.strip().split(",")]

    x = pd.DataFrame(ad.X.todense())
    x.columns = [x.split("-")[0] for x in ad.var.index]
    x.index = ad.obs.index

    # we need to remove cells with zero hashing counts
    # https://github.com/klarman-cell-observatory/demuxEM/pull/8/commits/fc5c4ebb28f01bce222774b5a035342d579b55bc
    x = x.loc[:, x.sum(axis=0)>0]

    outpath = outfile.replace(".sentinel",".gz")

    x.to_csv(outpath, index=True, index_label="Antibody")


    IOTools.touch_file(outfile)


@transform(hashCountCSV,
           regex(r'.*/.*/(.*).demuxEM.hash.count.csv.sentinel'),
           r"dehash.dir/demuxEM.dir/\1.demuxEM.sentinel")
def demuxEM(infile, outfile):
    '''
    Run demuxEM
    '''

    t = T.setup(infile, outfile, PARAMS, make_outdir = False)
    
    library_id = os.path.basename(infile).replace(".demuxEM.hash.count.csv.sentinel","")

    h5_file = os.path.join("api/counts/unfiltered",
                       library_id,
                       "h5/data.h5")

    hash_count_file = infile.replace(".sentinel",".gz")

    statement = '''demuxEM --generate-diagnostic-plots
                           --min-num-genes=100
                           --min-num-umis=100
                           --min-signal-hashtag=10
                            %(h5_file)s
                            %(hash_count_file)s
                            dehash.dir/demuxEM.dir/%(library_id)s
                 '''

    P.run(statement, **t.resources)

    IOTools.touch_file(outfile)


@transform(demuxEM,
           regex(r'.*/.*/(.*).demuxEM.sentinel'),
           r"dehash.dir/demuxEM.dir/\1.parse.demuxEM.sentinel")
def parseDemuxEM(infile, outfile):
        '''
        Parse the ridiculous output format from demuxEM.
        '''

        import pegasusio as pegio

        demuxEM_dir = os.path.dirname(infile)
        library_id = os.path.basename(infile).replace(".demuxEM.sentinel","")

        x = pegio.read_input(os.path.join(demuxEM_dir,
                                             library_id  + '.out.demuxEM.zarr.zip'))

        outpath = outfile.replace(".parse.demuxEM.sentinel",".tsv.gz")

        # fix the barcode.
        rp = pd.DataFrame(x.obsm["raw_probs"])

        rp.columns = [ "raw_prob_" + str(x) for x in rp.columns ]

        rp.index = x.obs.index

        out = pd.concat([x.obs[["demux_type","assignment"]], rp], axis=1)

        out.index = [ x + "-1" for x in out.index ]

        out.columns = [ "demuxEM_" + x for x in out.columns ]
        
        out["library_id"] = library_id 

        out.to_csv(outpath, index_label="barcode", sep="\t")

        IOTools.touch_file(outfile)


@merge(parseDemuxEM,
       "dehash.dir/demuxEM.api.sentinel")
def demuxemAPI(infiles, outfile):
    '''
    Register the demuxEM results on the API
    '''

    file_set={}

    for x in infiles:

        library_id = os.path.basename(x)[:-len(".parse.demuxEM.sentinel")]

        tsv_path = os.path.join(os.path.dirname(x),
                                library_id + ".tsv.gz")

        file_set[library_id] = {"path": tsv_path,
                                "description":"Parsed demuxEM results for " +\
                                library_id,
                                "format":"tsv"}

    x = T.api("dehash")

    x.define_dataset(analysis_name="demuxEM",
              data_subset="filtered",
              file_set=file_set,
              analysis_description="per library tables of demuxEM results")

    x.register_dataset()


# ########################################################################### #
# ##################### full target: to run all tasks ####################### #
# ########################################################################### #

@follows(gmmAPI)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
