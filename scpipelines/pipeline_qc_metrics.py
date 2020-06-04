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
      out_sentinel = "/".join(["qc_metrics", out_sample])
      infile = None
      yield(infile, out_sentinel)

@follows(mkdir("qc_metrics"))
@files(qc_metrics_jobs)
def calculate_qc_metrics(infile, outfile):
    ''' '''
  
    # Get cellranger directory
    sample_name = outfile.split("/")[-1].replace("_qcmetrics.tsv.gz", "")
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)
    cellranger_dir = samples.loc[sample_name, "path"]

    # Get genesets file
    if PARAMS["calculate_qc_metrics_geneset_file"] == "none" or PARAMS["calculate_qc_metrics_geneset_file"] == None:
      genesets_file = ""
    else:
      genesets_file = PARAMS["calculate_qc_metrics_geneset_file"]
      genesets_file = '''--genesets_file=%(genesets_file)s''' % locals()
 
    # Other settings 
    job_memory = PARAMS["resources_memory_high"]
    job_threads = PARAMS["resources_cores"]
    log_file = outfile.replace(".tsv.gz", ".log")
  
    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/qc_calculate_qc_metrics.R
                 --cellranger_dir=%(cellranger_dir)s
                 --numcores=%(job_threads)s
                 --log_filename=%(log_file)s
                 --outfile=%(outfile)s
                 %(genesets_file)s
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
      out_sentinel = "/".join(["scrublet", out_sample])
      infile = None
      yield(infile, out_sentinel)

@follows(mkdir("scrublet"))
@files(qc_doublet_scoring_jobs)
def run_scrublet(infile, outfile):
    ''' '''
  
    # Get cellranger directory
    sample_name = outfile.split("/")[-1].replace("_scrublet.sentinel", "")
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)
    cellranger_dir = samples.loc[sample_name, "path"]

    # Scrublet parameters
    expected_doublet_rate = PARAMS["scrublet_expected_doublet_rate"]
    min_counts = PARAMS["scrublet_min_counts"]
    min_cells = PARAMS["scrublet_min_cells"]
    min_gene_variability_pctl = PARAMS["scrublet_min_gene_variability_pctl"]
    n_prin_comps = PARAMS["scrublet_n_prin_comps"]
 
    # Other settings 
    job_memory = PARAMS["resources_memory_high"]
    job_threads = PARAMS["resources_cores"]
    log_file = outfile.replace(".sentinel", ".log")
    outdir = Path(outfile).parent
  
    # Formulate and run statement
    statement = '''python %(code_dir)s/python/run_scrublet.py
                   --cellranger_dir=%(cellranger_dir)s
                   --id=%(sample_name)s
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

@follows(calculate_qc_metrics, run_scrublet)
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
