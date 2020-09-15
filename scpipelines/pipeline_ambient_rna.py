"""==================
Pipeline ambient rna
=====================

Overview
========
This pipeline performs the following steps:
* Analyse the ambient RNA profile in each input (eg. channel's or sample's raw cellrange matrices)
* Compare ambient RNA profiles across inputs


Configuration
------------
The pipeline requires a configured :file:`pipeline.yml` file.
Default configuration files can be generated by executing:
   python <srcdir>/pipeline_ambient_rna.py config


Input files
-----------
An tsv file called 'input_samples.tsv' is required.
This file must have column names as explained below.
Must not include row names.
Add as many rows as iput channels/samples for analysis.

This file must have the following columns:

* sample_id - name used throughout. This could be the channel_pool id eg. A1
* path - path to the raw_matrix folder from cellranger count
* exp_batch - might or might not be useful. If not used, fill with "1"
* channel_id - might or might not be useful. If not used, fill with "1"
* seq_batch - might or might not be useful. If not used, fill with "1"
* (optional) blacklist - path to a file with cell_ids to blacklist

You can add any other columns as required, for example pool_id

Dependencies
------------
This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* R dependencies required in the r scripts

Pipeline output
===============
The pipeline returns:
* per-input html report and tables saved in a 'profile_per_input' folder
* ambient rna comparison across inputs saved in a 'profile_compare' folder

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
import yaml

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
        if not os.path.exists(p):
          raise ValueError('Input folder from cellranger run (outs/)'
                             ' does not exist.')
    IOTools.touch_file(outfile)

# ########################################################################### #
# ########################### Ambient RNA analysis ########################## #
# ########################################################################### #

# ------------------------------------------------------------------------
# Create output folder for each input (e.g channel, sample) in "per_input"

@follows(checkInputs)
def genClusterJobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    infile = None

    for sample in samples["sample_id"]:
        outfolder = "profile_per_input.dir/" + sample
        outfile = os.path.join(outfolder, "prep.sentinel")
        yield(infile, outfile)

@follows(checkInputs)
@files(genClusterJobs)
def prepFolders(infile, outfile):
    ''' Prepare folder structure for samples '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    IOTools.touch_file(outfile)

# --------------------------------------------------------
# Run ambient rna analysis per input (e.g channel, sample)

@transform(prepFolders,
           regex(r"profile_per_input.dir/(.*)/prep.sentinel"),
           r"profile_per_input.dir/\1/ambient_rna.sentinel")
def ambient_rna_per_input(infile, outfile):
    '''Explore count and gene expression profiles of ambient RNA droplets per input
    - The output is saved in profile_per_input.dir/<input_id>
    - The output consists on a html report and a ambient_genes.txt.gz file
    - See more details of the output in the ambient_rna_per_sample.R
    '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create options dictionary
    options = {}
    options["umi"] = int(PARAMS["ambientRNA_umi"])
    print(infile)
    sample_name = infile.split("/")[1].replace(".sample.dir", "")
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    samples.set_index("sample_id", inplace=True)
    options["cellranger_dir"] = samples.loc[sample_name ,"path"]
    options["outdir"] = outdir
    options["sample_name"] = sample_name

    # remove blacklisted cells if required
    if 'blacklist' in samples.columns:
        options["blacklist"] = samples.loc[sample_name, "blacklist"]

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir, "ambient_rna.yml"))
    print(task_yaml_file)
    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)
    output_dir = os.path.abspath(outdir)
    knit_root_dir = os.getcwd()
    fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")
    job_threads = PARAMS["resources_threads"]

    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/ambient_rna_per_sample.R
                   --task_yml=%(task_yaml_file)s
                   --log_filename=%(log_file)s
                '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)


# ------------------------------------------------------------------
# Compare ambient rna profiles from all inputs (e.g channel, sample)

@merge(ambient_rna_per_input,
       "profile_compare.dir/ambient_rna_compare.sentinel")
def ambient_rna_compare(infile, outfile):
    '''Compare the expression of top ambient RNA genes across inputs
    - The output is saved in profile_compare.dir
    - Output includes and html report and a ambient_rna_profile.tsv
    - See more details of the output in the ambient_rna_compare.R
    '''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create options dictionary
    options = {}
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    sample_id = samples.sample_id.tolist()
    sample_indir = [ "profile_per_input.dir/" + s for s in sample_id ]
    sample_indir = ",".join(sample_indir)
    sample_id = ",".join(sample_id)
    options["sample_indir"] = sample_indir
    options["sample_id"] = sample_id
    options["sample_table"] = "input_samples.tsv"
    options["outdir"] = outdir

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir, "ambient_rna_compare.yml"))
    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)
    output_dir = os.path.abspath(outdir)
    knit_root_dir = os.getcwd()
    fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")
    job_threads = PARAMS["resources_threads"]

    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/ambient_rna_compare.R
                   --task_yml=%(task_yaml_file)s
                   --log_filename=%(log_file)s
                '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)

# ----------------------
# Generic pipeline tasks

@follows(ambient_rna_compare)
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
