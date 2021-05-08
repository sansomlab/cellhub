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

def per_input(infile, input_libraries, outfile, PARAMS):
    # Create options dictionary

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    options = {}
    options["umi"] = int(PARAMS["ambientRNA_umi"])

    library_name = infile.split("/")[2].replace(".library.dir", "")
    libraries = pd.read_csv(input_libraries, sep='\t')
    libraries.set_index("library_id", inplace=True)
    options["cellranger_dir"] = libraries.loc[library_name ,"raw_path"]
    options["outdir"] = outdir
    options["library_name"] = library_name

    # remove blacklisted cells if required
    if 'blacklist' in libraries.columns:
        options["blacklist"] = libraries.loc[library_name, "blacklist"]

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
    statement = '''Rscript %(code_dir)s/R/ambient_rna_per_library.R
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
    libraries = pd.read_csv(PARAMS["input_libraries"], sep='\t')
    library_id = libraries.library_id.tolist()
    library_indir = [ "ambient.rna.dir/profile_per_input.dir/" + s for s in library_id ]
    library_indir = ",".join(library_indir)
    library_id = ",".join(library_id)
    options["library_indir"] = library_indir
    options["library_id"] = library_id
    options["library_table"] = "input_libraries.tsv"
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
