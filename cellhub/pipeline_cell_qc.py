'''
===================
pipeline_cell_qc.py
===================

Overview
========
This pipeline performs the following steps:

* Calculates per-cell QC metrics: ngenes, total_UMI, pct_mitochondrial, pct_ribosomal, pct_immunoglobin, pct_hemoglobin, and any specified geneset percentage
* Runs scrublet to calculate per-cell doublet score


Configuration
-------------
The pipeline requires a configured :file:`pipeline.yml` file. Default configuration files can be generated by executing: ::

   python <srcdir>/pipeline_cell_qc.py config


Input files
-----------
A tsv file called 'libraries.tsv' is required.
This file must have column names as explained below.
Must not include row names.
Add as many rows as input channels/librarys for analysis.
This file must have the following columns:
* library_id - name used throughout. This could be the channel_pool id eg. A1
* path - path to the filtered_matrix folder from cellranger count


Dependencies
------------
This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* R dependencies required in the r scripts


Pipeline output
---------------
The pipeline returns:
* qcmetrics.dir folder with per-input qcmetrics.tsv.gz table
* scrublet.dir folder with per-input scrublet.tsv.gz table

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

import cellhub.tasks as T


# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = T.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(T.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]


# ----------------------------- Preflight checks ----------------------------- #

if len(sys.argv) > 1:
    if sys.argv[1] == "make":
        inputs = glob.glob("api/counts/filtered/*/mtx/matrix.mtx.gz")
        if len(inputs) == 0:
            raise ValueError('No input files found on api/counts. (Counts from '
                             'the upstream pipeline can be registered with '
                            '"cellhub [upstream_pipeline_name] make useCounts")')


# ------------------------------ Pipeline Tasks ------------------------------ #

@follows(mkdir("cell.qc.dir"))
@transform(glob.glob("api/counts/filtered/*/mtx/matrix.mtx.gz"),
           regex(r".*/.*/.*/(.*)/mtx/matrix.mtx.gz"),
           r"cell.qc.dir/qcmetric.dir/\1.sentinel")
def qcmetrics(infile, outfile):
    '''This task will run R/calculate_qc_metrics.R,
    It uses the input_libraries.tsv to read the path to the cellranger directory for each input
    Ouput: creates a cell.qc.dir folder and a library_qcmetrics.tsv.gz table per library/channel
    For additional input files check the calculate_qc_metrics pipeline.yml sections:
    - Calculate the percentage of UMIs for genesets provided
    - Label barcodes as True/False based on whether they are part or not of a set of lists of barcodes provided
    '''
    
    t = T.setup(infile, outfile, PARAMS, 
                memory=PARAMS["resources_job_memory"], 
                cpu=PARAMS["resources_threads"])

    # Get cellranger directory and id
    library_name = os.path.basename(outfile)[:-len(".sentinel")]
    cellranger_dir = os.path.dirname(infile)

    # Get genesets file
    if PARAMS["calculate_qc_metrics_geneset_file"] == "none" or PARAMS["calculate_qc_metrics_geneset_file"] == None:
      genesets_file = ""
    else:
      genesets_file = PARAMS["calculate_qc_metrics_geneset_file"]
      genesets_file = '''--genesets_file=%(genesets_file)s''' % locals()

    # Get file with files having barcodes to label as 'True' in output dataframe
    if PARAMS["calculate_qc_metrics_barcodes_to_label_as_True"] == "none" or PARAMS["calculate_qc_metrics_barcodes_to_label_as_True"] == None:
      barcodes_to_label_as_True = ""
    else:
      barcodes_to_label_as_True = PARAMS["calculate_qc_metrics_barcodes_to_label_as_True"]
      barcodes_to_label_as_True = '''--barcodes_to_label_as_True=%(barcodes_to_label_as_True)s''' % locals()


    # Formulate and run statement
    statement = '''Rscript %(cellhub_code_dir)s/R/scripts/qc_metrics.R
                 --cellranger_dir=%(cellranger_dir)s
                 --library_id=%(library_name)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(out_file)s.tsv.gz
                 %(genesets_file)s
                 %(barcodes_to_label_as_True)s
              ''' % dict(PARAMS, **t.var, **locals())
              
    P.run(statement, **t.resources)

    # Create sentinel file
    IOTools.touch_file(outfile)


@merge(qcmetrics,
       "cell.qc.dir/qcmetrics.dir/api.sentinel")
def qcmetricsAPI(infiles, outfile):
    '''
    Add the QC metrics results to the API
    '''

    file_set={}

    for libqc in infiles:

        tsv_path = libqc.replace(".sentinel",".tsv.gz")
        library_id = os.path.basename(tsv_path)

        file_set[library_id] = {"path": tsv_path,
                                "description":"qcmetric table for library " +\
                                library_id,
                                "format":"tsv"}

    x = T.api("cell.qc")

    x.define_dataset(analysis_name="qcmetrics",
              data_subset="filtered",
              file_set=file_set,
              analysis_description="per library tables of cell GEX qc statistics")

    x.register_dataset()



@follows(mkdir("cell.qc.dir"))
@transform(glob.glob("api/counts/filtered/*/mtx/matrix.mtx.gz"),
           regex(r".*/.*/.*/(.*)/mtx/matrix.mtx.gz"),
           r"cell.qc.dir/scrublet.dir/\1.sentinel")
def scrublet(infile, outfile):
    '''This task will run python/run_scrublet.py,
    It uses the input_libraries.tsv to read the path to the cellranger directory for each input
    Ouput: creates a scrublet.dir folder and a library_scrublet.tsv.gz table per library/channel
    It also creates a doublet score histogram and a double score umap for each library/channel
    Check the scrublet section in the pipeline.yml to specify other parameters
    '''

    t = T.setup(infile, outfile, PARAMS, memory="50G", cpu=3)

    library_name = os.path.basename(outfile)[:-len(".sentinel")]
    cellranger_dir = os.path.dirname(infile)

    if PARAMS["scrublet_subset"]:
        whitelist = libraries.loc[library_name, "whitelist"]
        subset_option = '''--keep_barcodes_file=%(whitelist)s''' %locals()
    else:
        subset_option = ''' '''

    # Scrublet parameters
    expected_doublet_rate = PARAMS["scrublet_expected_doublet_rate"]
    min_counts = PARAMS["scrublet_min_counts"]
    min_cells = PARAMS["scrublet_min_cells"]
    min_gene_variability_pctl = PARAMS["scrublet_min_gene_variability_pctl"]
    n_prin_comps = PARAMS["scrublet_n_prin_comps"]

    # Formulate and run statement
    statement = '''python %(cellhub_code_dir)s/python/qc_scrublet.py
                   --cellranger_dir=%(cellranger_dir)s
                   %(subset_option)s
                   --library_id=%(library_name)s
                   --expected_doublet_rate=%(expected_doublet_rate)s
                   --min_counts=%(min_counts)s
                   --min_cells=%(min_cells)s
                   --min_gene_variability_pctl=%(min_gene_variability_pctl)s
                   --n_prin_comps=%(n_prin_comps)s
                   --outdir=%(outdir)s
                   &> %(log_file)s
                ''' % dict(PARAMS, **t.var, **locals())

    P.run(statement, **t.resources)

    # Create sentinel file
    IOTools.touch_file(outfile)


@merge(scrublet,
       "cell.qc.dir/scrublet.dir/api.sentinel")
def scrubletAPI(infiles, outfile):
    '''
    Add the scrublet results to the API
    '''

    file_set={}

    for lib in infiles:

        tsv_path = lib.replace(".sentinel",".tsv.gz")
        library_id = os.path.basename(tsv_path)

        file_set[library_id] = {"path": tsv_path,
                                "description":"scrublet table for library " +\
                                library_id,
                                "format":"tsv"}

    x = T.api("cell.qc")

    x.define_dataset(analysis_name="scrublet",
              data_subset="filtered",
              file_set=file_set,
              analysis_description="per library tables of cell scrublet scores")

    x.register_dataset()


# ---------------------------------------------------
# Generic pipeline tasks

@follows(mkdir("cell.qc.dir"))
@files(None, "cell.qc.dir/plot.sentinel")
def plot(infile, outfile):
    '''Draw the pipeline flowchart'''

    pipeline_printout_graph ( "cell.qc.dir/pipeline_flowchart.svg",
                          "svg",
                          [full],
                          no_key_legend=True)

    pipeline_printout_graph ( "cell.qc.dir/pipeline_flowchart.png",
                          "png",
                          [full],
                          no_key_legend=True)

    IOTools.touch_file(outfile)


@follows(qcmetricsAPI, scrubletAPI, plot)
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
