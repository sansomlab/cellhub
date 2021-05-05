'''
Tasks for pipeline_ambient_rna.py
=================================

'''
import sys
import os
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools
from pathlib import Path
import pandas as pd
import yaml
from . import TASK

def per_input(infile, input_samples, outfile, PARAMS):
    # Create options dictionary

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    options = {}
    options["umi"] = int(PARAMS["ambientRNA_umi"])

    sample_name = infile.split("/")[2].replace(".sample.dir", "")
    samples = pd.read_csv(input_samples, sep='\t')
    samples.set_index("sample_id", inplace=True)
    options["cellranger_dir"] = samples.loc[sample_name ,"raw_path"]
    options["outdir"] = outdir
    options["sample_name"] = sample_name

    # remove blacklisted cells if required
    if 'blacklist' in samples.columns:
        options["blacklist"] = samples.loc[sample_name, "blacklist"]

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir, "ambient_rna.yml"))

    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)
    output_dir = os.path.abspath(outdir)
    knit_root_dir = os.getcwd()
    fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")
    job_threads = PARAMS["resources_threads"]

    job_threads, job_memory, r_memory = TASK.get_resources(
        memory=PARAMS["resources_job_memory"])

    # Formulate and run statement
    statement = '''Rscript %(code_dir)s/R/ambient_rna_per_sample.R
                   --task_yml=%(task_yaml_file)s
                   --log_filename=%(log_file)s
                '''

    P.run(statement)


def compare(infile, outfile, PARAMS):

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create options dictionary
    options = {}
    samples = pd.read_csv("input_samples.tsv", sep='\t')
    sample_id = samples.sample_id.tolist()
    sample_indir = [ "ambient.rna.dir/profile_per_input.dir/" + s for s in sample_id ]
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
