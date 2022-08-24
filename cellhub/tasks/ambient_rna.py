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

def per_input(infile, outfile, PARAMS):


    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    mtx_dir = os.path.dirname(infile)

    library_id = Path(outfile).parents[0]

    # Create options dictionary
    options = {}
    options["umi"] = int(PARAMS["ambientRNA_umi"])
    options["cellranger_dir"] = mtx_dir
    options["outdir"] = outdir
    options["library_name"] = library_id

    # remove blacklisted cells if required
    if PARAMS["blacklist"] is not None:
        options["blacklist"] = PARAMS["blacklist"]

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir,
                                                  "ambient_rna.yml"))

    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)

    # output_dir = os.path.abspath(outdir)
    # knit_root_dir = os.getcwd()
    # fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")

    # job_threads = PARAMS["resources_threads"]

    job_threads, job_memory, r_memory = TASK.get_resources(
        memory=PARAMS["resources_job_memory"], cpu=PARAMS["resources_threads"])

    # Formulate and run statement
    statement = '''Rscript %(cellhub_code_dir)s/R/scripts/ambient_rna_per_library.R
                   --task_yml=%(task_yaml_file)s
                   --log_filename=%(log_file)s
                '''

    P.run(statement)


def compare(infiles, outfile, PARAMS):

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    library_indir = ",".join([os.path.dirname(x) for x in infiles])
    library_id = ",".join([str(os.path.basename(Path(x).parents[0]))
                           for x in infiles])

    # Create options dictionary
    options = {}
    options["library_indir"] = library_indir
    options["library_id"] = library_id
    options["library_table"] = "input_libraries.tsv"
    options["outdir"] = outdir

    # Write yml file
    task_yaml_file = os.path.abspath(os.path.join(outdir, "ambient_rna_compare.yml"))
    with open(task_yaml_file, 'w') as yaml_file:
        yaml.dump(options, yaml_file)

    # output_dir = os.path.abspath(outdir)
    # knit_root_dir = os.getcwd()
    # fig_path =  os.path.join(output_dir, "fig.dir/")

    # Other settings
    log_file = outfile.replace("sentinel","log")
    job_threads = PARAMS["resources_threads"]

    if ("G" in PARAMS["resources_job_memory"] or
        "M" in PARAMS["resources_job_memory"] ):
        job_memory = PARAMS["resources_job_memory"]

    # Formulate and run statement
    statement = '''Rscript %(cellhub_code_dir)s/R/scripts/ambient_rna_compare.R
                   --task_yml=%(task_yaml_file)s
                   --log_filename=%(log_file)s
                '''
    P.run(statement)
