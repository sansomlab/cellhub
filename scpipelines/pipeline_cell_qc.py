"""==================
Pipeline Cell QC
=====================


Overview
========
This pipeline performs the following steps:
* Calculates per-cell QC metrics: ngenes, total_UMI, pct_mitochondrial,
  pct_ribosomal, pct_immunoglobin, pct_hemoglobin, and any specified geneset percentage
* Runs scrublet to calculate per-cell doublet score


Configuration
------------
The pipeline requires a configured :file:`pipeline.yml` file.
Default configuration files can be generated by executing:
   python <srcdir>/pipeline_cell_qc.py config


Input files
-----------
A tsv file called 'input_samples.tsv' is required.
This file must have column names as explained below.
Must not include row names.
Add as many rows as input channels/samples for analysis.
This file must have the following columns:
* sample_id - name used throughout. This could be the channel_pool id eg. A1
* path - path to the filtered_matrix folder from cellranger count


Dependencies
------------
This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* R dependencies required in the r scripts


Pipeline output
===============
The pipeline returns:
* qcmetrics.dir folder with per-input qcmetrics.tsv.gz table
* scrublet.dir folder with per-input scrublet.tsv.gz table

Code
====
"""
from ruffus import *
from ruffus.combinatorics import *
import sys
import os
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools
from pathlib import Path
import pandas as pd

# -------------------------- < parse parameters > --------------------------- #

# Load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# Set the location of the cellhub code directory
if "code_dir" not in PARAMS.keys():
    PARAMS["code_dir"] = Path(__file__).parents[1]
else:
    raise ValueError("Could not set the location of the code directory")

# ----------------------- < pipeline configuration > ------------------------ #

# handle pipeline configuration
if len(sys.argv) > 1:
        if(sys.argv[1] == "config") and __name__ == "__main__":
                    sys.exit(P.main(sys.argv))


# ########################################################################### #
# ######## Check input samples file and that the input exists ############### #
# ########################################################################### #

@originate("input.check.sentinel")
def checkInputs(outfile):
    '''Check that input_samples.tsv exists and the path given in the file
       is a valid directorys. '''

    if not os.path.exists("input_samples.tsv"):
        raise ValueError('File specifying the input samples is not present.'
                         'The file needs to be named "input_samples.tsv" ')

    samples = pd.read_csv("input_samples.tsv", sep='\t')
    for p in samples["path"]:
        print(p)
        if not os.path.exists(p):
          raise ValueError('Input folder from cellranger run (outs/)'
                             ' does not exist.')
    IOTools.touch_file(outfile)

# ############################################# #
# ######## Calculate QC metrics ############### #
# ############################################# #

@follows(checkInputs)
def qc_metrics_jobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)

    for sample_name in samples.index:
      out_sample = "_".join([sample_name, "qcmetrics.tsv.gz"])
      out_sentinel = "/".join(["qcmetrics.dir", out_sample])
      infile = None
      yield(infile, out_sentinel)

@follows(mkdir("qcmetrics.dir"))
@files(qc_metrics_jobs)
def calculate_qc_metrics(infile, outfile):
    '''This task will run R/calculate_qc_metrics.R,
    It uses the input_samples.tsv to read the path to the cellranger directory for each input
    Ouput: creates a qcmetrics.dir folder and a sample_qcmetrics.tsv.gz table per sample/channel
    For additional input files check the calculate_qc_metrics pipeline.yml sections:
    - Calculate the percentage of UMIs for genesets provided
    - Label barcodes as True/False based on whether they are part or not of a set of lists of barcodes provided
    '''

    # Get cellranger directory and id
    sample_name = outfile.split("/")[-1].replace("_qcmetrics.tsv.gz", "")
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)
    cellranger_dir = samples.loc[sample_name, "filt_path"]

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

    # Other settings
    job_threads = PARAMS["resources_threads"]
    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    log_file = outfile.replace(".tsv.gz", ".log")

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/calculate_qc_metrics.R
                 --cellranger_dir=%(cellranger_dir)s
                 --sample=%(sample_name)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(outfile)s
                 %(genesets_file)s
                 %(barcodes_to_label_as_True)s
              '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)

@follows(checkInputs)
def qc_reports_jobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)

    for sample_name in samples.index:
      out_sample = "_".join([sample_name, "qcmetrics_report.pdf"])
      out_sentinel = "/".join(["qcmetrics.dir/reports", out_sample])
      infile = None
      yield(infile, out_sentinel)

@follows(mkdir("qcmetrics.dir/reports"))
@files(qc_reports_jobs)
def build_qc_reports(infile, outfile):
    '''This task will run R/build_qc_mapping_report.R,
    It expects three files in the input directory barcodes.tsv.gz, 
    features.tsv.gz, and matrix.mtx.gz
    Ouput: creates a sample_qcmetrics_report.pdf table per input folder
    '''
    # Get cellranger directory and id
    sample_name = outfile.split("/")[-1].replace("_qcmetrics.dir/reports", "")
    sample_name = sample_name.replace("_qcmetrics_report.pdf", "")

    # cellranger filtered output
    cellranger_dir = "-".join([sample_name, "count/outs/filtered_feature_bc_matrix"])

    # Other settings
    job_threads = PARAMS["resources_threads"]
    job_threads = PARAMS["resources_threads"]
    job_memory = "30G"

    log_file = outfile.replace(".pdf", ".log")

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/build_qc_mapping_reports.R 
                --tenxfolder=%(cellranger_dir)s 
                --sample_id=%(sample_name)s
                --specie="hg"
                --outfolder="qcmetrics.dir/reports"
                &> %(log_file)s
              '''
    P.run(statement)





# ############################################# #
# ######## Calculate doublet scores ########### #
# ############################################# #

@follows(checkInputs)
def qc_doublet_scoring_jobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)

    for sample_name in samples.index:
      out_sample = "_".join([sample_name, "scrublet.sentinel"])
      out_sentinel = "/".join(["scrublet.dir", out_sample])
      infile = None
      yield(infile, out_sentinel)

@follows(mkdir("scrublet.dir"))
@files(qc_doublet_scoring_jobs)
def run_scrublet(infile, outfile):
    '''This task will run python/run_scrublet.py,
    It uses the input_samples.tsv to read the path to the cellranger directory for each input
    Ouput: creates a scrublet.dir folder and a sample_scrublet.tsv.gz table per sample/channel
    It also creates a doublet score histogram and a double score umap for each sample/channel
    Check the scrublet section in the pipeline.yml to specify other parameters
    '''

    # Get cellranger directory
    sample_name = outfile.split("/")[-1].replace("_scrublet.sentinel", "")
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)
    cellranger_dir = samples.loc[sample_name, "filt_path"]

    if PARAMS["scrublet_subset"]:
        whitelist = samples.loc[sample_name, "whitelist"]
        subset_option = '''--keep_barcodes_file=%(whitelist)s''' %locals()
    else:
        subset_option = ''' '''

    # Scrublet parameters
    expected_doublet_rate = PARAMS["scrublet_expected_doublet_rate"]
    min_counts = PARAMS["scrublet_min_counts"]
    min_cells = PARAMS["scrublet_min_cells"]
    min_gene_variability_pctl = PARAMS["scrublet_min_gene_variability_pctl"]
    n_prin_comps = PARAMS["scrublet_n_prin_comps"]

    # Other settings
    job_threads = PARAMS["resources_threads"]
    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]
    
    job_threads=3
    job_memory="50G"
    log_file = outfile.replace(".sentinel", ".log")
    outdir = Path(outfile).parent

    # Formulate and run statement
    statement = '''python %(code_dir)s/python/run_scrublet.py
                   --cellranger_dir=%(cellranger_dir)s
                   %(subset_option)s
                   --sample=%(sample_name)s
                   --expected_doublet_rate=%(expected_doublet_rate)s
                   --min_counts=%(min_counts)s
                   --min_cells=%(min_cells)s
                   --min_gene_variability_pctl=%(min_gene_variability_pctl)s
                   --n_prin_comps=%(n_prin_comps)s
                   --outdir=%(outdir)s
                   &> %(log_file)s
                '''

    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)


# ---------------------------------------------------
# Generic pipeline tasks

@follows(calculate_qc_metrics, build_qc_reports, run_scrublet)
def full():
    '''
    Run the full pipeline.
    '''
    pass


pipeline_printout_graph ( "pipeline_flowchart.svg",
                          "svg",
                          [full],
                          no_key_legend=True)

pipeline_printout_graph ( "pipeline_flowchart.png",
                          "png",
                          [full],
                          no_key_legend=True)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
