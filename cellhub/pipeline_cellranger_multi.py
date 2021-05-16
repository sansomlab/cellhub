'''
=========================
Pipeline Cellranger Multi
=========================


Overview
========

This pipeline performs the following functions:

* Alignment and quantitation (using cellranger count or cellranger multi)

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_cellranger_multi.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_cellranger_multi.py config


Dependencies
------------

This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* cellranger: https://support.10xgenomics.com/single-cell-gene-expression/


Pipeline output
===============

The pipeline returns:
* the output of cellranger multi

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
from cellhub.tasks import templates
from cellhub.tasks import resources
from cellhub.tasks import TASK

import cellhub.tasks.control as C
import cellhub.tasks.cellranger as cellranger
import cellhub.tasks.api as api

# Override function to collect config files
P.control.write_config_files = C.write_config_files

# -------------------------- < parse parameters > --------------------------- #

# load options from the yml file
parameter_file = C.get_parameter_file(__file__, __name__)
PARAMS = P.get_parameters(parameter_file)

# set the location of the pipeline code directory
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


# ----------------------- < helper functions > ------------------------ #


@files(None, "task.summary.table.tex")
def taskSummary(infile, outfile):
    '''Make a summary of optional tasks that will be run'''

    tasks, run = [], []

    for k,v in PARAMS.items():
        if k.startswith("run_"):
            tasks.append(k[4:])
            run.append(str(v))

    tab = pd.DataFrame(list(zip(tasks,run)),columns=["task","run"])

    tab.to_latex(buf=outfile, index=False)



# ########################################################################### #
# ################ Read parameters and create config file ################### #
# ########################################################################### #

@active_if(PARAMS["input"] == "mkfastq")
@follows(mkdir("cellranger.multi.dir"))
@originate("cellranger.multi.dir/config.sentinel")
def makeConfig(outfile):
    '''Read parameters from yml file for the whole experiment and save config files as csv.'''

    # check if references exist
    if PARAMS["run_gene-expression"]:
        gexref = PARAMS["gene-expression_reference"]
        if gexref is None:
            raise ValueError('"gene-expression_reference" parameter not set'
                             ' in file "pipeline.yml"')

        if not os.path.exists(gexref):
            raise ValueError('The specified "gene-expression_reference"'
                             ' file does not exist')
    else:
        pass

    if PARAMS["run_feature"]:
        featureref = PARAMS["feature_reference"]
        if featureref is None:
            raise ValueError('"feature_reference" parameter not set'
                             ' in file "pipeline.yml"')

        if not os.path.exists(featureref):
            raise ValueError('The specified "feature_reference"'
                             ' file does not exist')
    else:
        pass

    if PARAMS["run_vdj"]:
        vdjref = PARAMS["vdj_reference"]
        if vdjref is None:
            raise ValueError('"vdj_reference" parameter not set'
                             ' in file "pipeline.yml"')

        if not os.path.exists(vdjref):
            raise ValueError('The specified "vdj_reference"'
                             ' file does not exist')
    else:
        pass


    # read parameters for gex
    section, param = [], []

    if PARAMS["run_gene-expression"]:
        for k,v in PARAMS.items():
            if k.startswith("gene-expression_"):
                if v is not None:
                    section.append(k[16:])
                    param.append(str(v))

        df_gex = pd.DataFrame(list(zip(section,param)),columns=["[gene-expression]",""])
    else:
        pass

    # read parameters for feature
    section, param = [], []

    if PARAMS["run_feature"]:
        for k,v in PARAMS.items():
            if k.startswith("feature_"):
                if v is not None:
                    section.append(k[8:])
                    param.append(str(v))

        df_feature = pd.DataFrame(list(zip(section,param)),columns=["[feature]",""])
    else:
        pass

    # read parameters for vdj
    section, param = [], []

    if PARAMS["run_vdj"]:
        for k,v in PARAMS.items():
            if k.startswith("vdj_"):
                if v is not None:
                    section.append(k[4:])
                    param.append(str(v))

        df_vdj = pd.DataFrame(list(zip(section,param)),columns=["[vdj]",""])
    else:
        pass

    # read parameters for libraries:
    lib_params = PARAMS["libraries"]
    library_ids = list(lib_params.keys())

    for library_id in library_ids:

        # Save subsections of parameters in config files specific for each sample
        # (data.dir/sample01.csv data.dir/sample02.csv etc)
        libsample_params = PARAMS["libraries_" + library_id]

        filename = "cellranger.multi.dir/" + library_id + ".csv"

        lib_df = pd.DataFrame(libsample_params)
        print(lib_df)

        lib_df.drop('description', axis=1, inplace=True)

        lib_columns = list(lib_df)

        smp_df = pd.DataFrame()
        for i in lib_columns:
            tmp = lib_df[i].str.split(',', expand=True)
            smp_df = smp_df.append(tmp.T)

            # filter out gex rows from libraries table if run_gene-expression = false
            mask = smp_df.feature_types == 'Gene Expression'
            if PARAMS["run_gene-expression"]:
                df_filt = smp_df
            else:
                df_filt = smp_df[~mask]

            # filter out feature rows from libraries table if run_feature = false
            mask = df_filt.feature_types == 'Antibody Capture'
            if PARAMS["run_feature"]:
                df_filt = df_filt
            else:
                df_filt = df_filt[~mask]

            # filter out vdj rows from libraries table if run_vdj = false
            mask = df_filt.feature_types == 'VDJ-B'
            if PARAMS["run_vdj"]:
                df_filt = df_filt
            else:
                df_filt = df_filt[~mask]


        # but I need to add different headers for each subsection, so I stream each table individually.
        with open(filename, 'a') as csv_stream:

            if PARAMS["run_gene-expression"]:
                csv_stream.write('[gene-expression]\n')
                df_gex.to_csv(csv_stream, header=False, index=False)
                csv_stream.write('\n')
            else:
                pass

            if PARAMS["run_feature"]:
                csv_stream.write('[feature]\n')
                df_feature.to_csv(csv_stream, header=False, index=False)
                csv_stream.write('\n')
            else:
                pass

            if PARAMS["run_vdj"]:
                csv_stream.write('[vdj]\n')
                df_vdj.to_csv(csv_stream, header=False, index=False)
                csv_stream.write('\n')
            else:
                pass

            csv_stream.write('[libraries]\n')
            df_filt.to_csv(csv_stream, header=True, index=False)
            csv_stream.write('\n')

    IOTools.touch_file(outfile)

# ########################################################################### #
# ############################ run cellranger multi ######################### #
# ########################################################################### #

@follows(makeConfig)
@transform("cellranger.multi.dir/*.csv",
           regex(r".*/([^.]*).*.csv"),
           r"cellranger.multi.dir/\1-cellranger.multi.sentinel")
def cellrangerMulti(infile, outfile):
    '''
    Execute the cellranger multi pipleline for first sample.
    '''

    # read id_tag from file name
    config_path = os.path.basename(infile)
    sample_basename = os.path.basename(infile)
    sample_name_sections = sample_basename.split(".")
    id_tag = sample_name_sections[0]


    #set the maximum number of jobs for cellranger
    max_jobs = PARAMS["cellranger_maxjobs"]

    ## send one job script to slurm queue which arranges cellranger run
    ## hard-coded to ensure enough resources
    job_threads = 6
    job_memory = "24G"

    log_file = id_tag + ".log"

    # this statement is to run in slurm mode
    statement = '''cd cellranger.multi.dir;
                    cellranger multi
	    	    --id %(id_tag)s
                    --csv=%(config_path)s
                    --jobmode=slurm
                    --maxjobs=%(max_jobs)s
		    --nopreflight
                    &> %(log_file)s
                 '''

    P.run(statement)
    IOTools.touch_file(outfile)

@transform(cellrangerMulti,
           regex(r"(.*)/(.*)-cellranger.multi.sentinel"),
           r"\1/out.dir/\2/post.process.matrices.sentinel")
def postProcessMatrices(infile, outfile):
    '''
    Post-process the cellranger multi matrices to split the
    counts for the GEX, ADT and HTO modalities into seperate
    market matrices.

    * A critical function of the pre-processing is to reformat
      the cellbarcodes to the "UMI-1-LIBRARY_ID" format.

    cellranger.multi.dir folder layout is

    (1) unfiltered outputs
    ----------------------
    library_id/outs/multi/count/raw_feature_bc_matrix/
    library_id/outs/multi/vdj_b/

    (2) filtered outputs
    --------------------
    library_id/outs/per_sample_outs/sample|library_id/count/sample_feature_bc_matrix
    library_id/outs/per_sample_outs/sample|library_id/vdj_b

    Notes
    -----
    - the feature_bc_matrix can contain GEX, ADT and HTO
    - vdj_b = BCR sequencing
    - vdj_t ?!? presumably this is what the TCR folder will look like but we have not
      had any datasets yet.

    Task overview
    -------------
    (i) split the raw counts into seperate GEX, HTO and ADT matrices
    (ii) link in the raw vdj
    (iii) for each sample, do (i) and (ii)

    Outputs
    -------

    post.processed.dir/unfiltered/gex/
    post.processed.dir/unfiltered/ADT/
    post.processed.dir/unfiltered/HTO/
    post.processed.dir/unfiltered/vdj_b/
    post.processed.dir/unfiltered/vdj_t/

    post.processed.dir/per_sample/sample|library_id/gex/
    post.processed.dir/per_sample/sample|library_id/ADT/
    post.processed.dir/per_sample/sample|library_id/HTO/
    post.processed.dir/per_sample/sample|library_id/vdj_b/
    post.processed.dir/per_sample/sample|library_id/vdj_t/
    '''

    out_dir = os.path.dirname(outfile)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    library_id = os.path.basename(infile).split("-cellranger.multi")[0]

    # 1. deal with unfiltered count data
    matrix_location = os.path.join("cellranger.multi.dir", library_id,
                                   "outs/multi/count/raw_feature_bc_matrix")

    output_location = os.path.join("cellranger.multi.dir/out.dir/",
                                   library_id, "unfiltered")

    cellranger.get_counts(matrix_location, output_location,
                          library_id)

    for vdj_type in ["b", "t"]:

        vdj_location = os.path.join("cellranger.multi.dir", library_id,
                                    "outs/multi/vdj_" + vdj_type)

        vdj_out_location = os.path.join("cellranger.multi.dir/out.dir/",
                                         library_id, "unfiltered/vdj" + vdj_type)

    if os.path.exists(vdj_location):
        os.symlink(vdj_location, vdj_b_out_location)


    # 2. deal with per sample libraries

    per_sample_loc = os.path.join("cellranger.multi.dir",
                                  library_id,
                                  "outs/per_sample_outs/")

    per_sample_dirs = glob.glob(per_sample_loc + "*")

    for per_sample_dir in per_sample_dirs:

        matrix_location = os.path.join(per_sample_dir,
                                       "count/sample_feature_bc_matrix")

        sample_id = os.path.basename(per_sample_dir)

        output_location = os.path.join("cellranger.multi.dir/out.dir/",
                                       library_id,
                                       sample_id)

        cellranger.get_counts(matrix_location, output_location,
                              library_id)

        for vdj_type in ["b", "t"]:

            vdj_location = os.path.join(per_sample_dir,
                                        "vdj_" + vdj_type)

            vdj_out_location = os.path.join("cellranger.multi.dir/out.dir/",
                                            library_id,
                                            sample_id,
                                            "vdj_" + vdj_type)

            if os.path.exists(vdj_location):
                os.symlink(vdj_location, vdj_b_out_location)


    IOTools.touch_file(outfile)



@transform(postProcessMatrices,
           regex(r"(.*)/out.dir/(.*)/post.process.matrices.sentinel"),
           r"\1/out.dir/\2/api.register.sentinel")
def API(infile, outfile):
    '''
    Register the outputs on the service endpoint
    '''

    # 1. register the GEX, ADT and HTO count matrices

    x = api.register("cellranger.multi")

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

    library_id = outfile.split("/")[-2]

    source_loc = os.path.dirname(infile)

    for data_subset in ["unfiltered", "filtered"]:


        for modality in ["GEX", "ADT", "HTO"]:

            if data_subset == "filtered":
                subset_dir = library_id
            else:
                subset_dir = data_subset

            mtx_loc = os.path.join(source_loc,
                                   subset_dir,
                                   modality)

            if os.path.exists(mtx_loc):

                mtx_x = mtx_template.copy()
                mtx_x["barcodes"]["path"] = os.path.join(mtx_loc, "barcodes.tsv.gz")
                mtx_x["features"]["path"] = os.path.join(mtx_loc, "features.tsv.gz")
                mtx_x["matrix"]["path"] =  os.path.join(mtx_loc, "matrix.mtx.gz")

                x.dataset(analysis_name=modality,
                          data_subset=data_subset,
                          data_id=library_id,
                          file_set=mtx_x,
                          analysis_description="unfiltered cellranger count GEX output",
                          file_format="10X-market-matrix")

                print("---------")
                x.report()
                x.deposit()





# @follows(makeConfig)
# @merge("cellranger.multi.dir/*.csv",
#         "cellranger.multi.dir/libraries.tsv")
# def makeLibraryTable(library_files, outfile):
#     # Build the path to the log file

#     library_names = []

#     for library_file in library_files:
#         library_name = os.path.basename(library_file)
#         library_names.append(library_name)

#     libraries = ','.join(library_names)

#     job_threads = 2
#     job_memory = "2000M"

#     log_file = outfile.replace(".tsv", ".log")
#     statement = '''Rscript %(code_dir)s/R/cellranger_library_table.R
#                     --outfile=%(outfile)s
#                     --librarydir=cellranger.multi.dir
#                     --libraryfiles=%(libraries)s
#                     &> %(log_file)s
#                 '''
#     P.run(statement)

#     IOTools.touch_file(outfile + ".sentinel")

#
# ---------------------------------------------------
# Generic pipeline tasks

@follows(cellrangerMulti, API)
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
