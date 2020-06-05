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

# ############################################# #
# ######## Ambient RNA ############### #
# ############################################# #

# Make sample folders first
@follows(checkInputs)
def genClusterJobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    infile = None

    for sample in samples["sample_id"]:
        outfolder = "per_input/" + sample
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


@transform(prepFolders,
           regex(r"per_input/(.*)/prep.sentinel"),
           r"per_input/\1/ambient_rna.sentinel")
def runAmbientRNA(infile, outfile):
    '''Explore count and gene expression profiles of ambient RNA droplets per sample'''

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
    job_memory = PARAMS["resources_memory_high"]

    # Formulate and run statement
    statement = '''Rscript -e "rmarkdown::render('%(code_dir)s/R/qc_ambient_rna_per_sample.R',
                   output_dir = '%(output_dir)s',
                   intermediates_dir = '%(output_dir)s',
                   knit_root_dir= '%(knit_root_dir)s',
                   params=list('task_yml' = '%(task_yaml_file)s',
                               'fig_path' = '%(fig_path)s',
                               'log_filename' = '%(log_file)s' ) )"
                '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)

@merge(runAmbientRNA,
       "summary.dir/ambient_rna_summary.sentinel")
def AmbientRNAsummarise(infile, outfile):
    '''Compare the expression of top ambient RNA genes across samples'''

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create options dictionary
    options = {}
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    sample_id = samples.sample_id.tolist()
    sample_indir = [ "per_input/" + s for s in sample_id ]
    sample_indir = ",".join(sample_indir)
    sample_id = ",".join(sample_id)
    options["sample_indir"] = sample_indir
    options["sample_id"] = sample_id
    options["sample_table"] = "input_samples.tsv"
    options["outdir"] = outdir
    if PARAMS["ambientRNA_plot_annotation"] != "none" :
      if PARAMS["ambientRNA_plot_annotation"] == "all" :
        cols = samples.columns.values.tolist()
        remove_cols = {'sample_id','path','blacklist', 'agg_id'}
        cols = [c for c in cols if c not in remove_cols]
        cols = ",".join(cols)
        options["plot_annotation"] = cols
      else:
        add_cols = PARAMS["ambientRNA_plot_annotation"]
        options["plot_annotation"] = "exp_batch,channel_id,seq_batch," + add_cols
    else:
      options["plot_annotation"] = "exp_batch,channel_id,seq_batch"

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir, "ambientRNA_summary.yml"))
    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)
    output_dir = os.path.abspath(outdir)
    knit_root_dir = os.getcwd()
    fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")
    job_memory = PARAMS["resources_memory_standard"]

    # Formulate and run statement
    statement = '''Rscript -e "rmarkdown::render('%(code_dir)s/R/qc_ambient_rna_all_samples.R',
                   output_dir = '%(output_dir)s',
                   intermediates_dir = '%(output_dir)s',
                   knit_root_dir= '%(knit_root_dir)s',
                   params=list('task_yml' = '%(task_yaml_file)s',
                               'fig_path' = '%(fig_path)s',
                               'log_filename' = '%(log_file)s' ) )"
                '''
    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)
    
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
