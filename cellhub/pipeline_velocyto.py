"""
=================
Pipeline Velocyto
=================

Overview
========

This pipeline performs the following steps:

* sort bam file by cell barcode
* estimate intronic and exonic reads using velocyto (on selected barcodes)

Usage
=====

See :doc:`Installation</Installation>` and :doc:`Usage</Usage>` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_velocity.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_velocity.py config


Input files
-----------

The pipeline is run from bam files generated by cellranger count.

The pipeline expects a tsv file containing the path to each cellranger
bam file (path) and the respective sample_id for each sample.
In addition a list of barcodes is required, this could be the filtered
barcodes from cellranger or a custom input (can be gzipped file).
Any further metadata can be added to the file. The required columns
are sample_id, barcodes and path.


Dependencies
------------

This pipeline requires:
* cgat-core: https://github.com/cgat-developers/cgat-core
* samtools
* veloctyo


Pipeline output
===============

The pipeline returns:
* a loom file with intronic and exonic reads for use in scvelo analysis

Code
====

"""
from ruffus import *
from pathlib import Path
import sys
import os
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools
import pandas as pd

import cellhub.tasks as T

# -------------------------- Pipeline Configuration -------------------------- #

# Override function to collect config files
P.control.write_config_files = T.write_config_files

# load options from the yml file
P.parameters.HAVE_INITIALIZED = False
PARAMS = P.get_parameters(T.get_parameter_file(__file__))

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]

# ------------------------------ Pipeline Tasks ------------------------------ #


@originate("input.check.sentinel")
def checkInputs(outfile):
    '''Check that input_samples.tsv exists and the path given in the file
       is a valid directorys. '''

    if not os.path.exists(PARAMS["input_samples"]):
        raise ValueError('File specifying the input samples is not present.'
                         'The file needs to be named "input_samples.txt" ')

    samples = pd.read_csv(PARAMS["input_samples"], sep='\t')
    for p in samples["path"]:
        print(p)
        if not os.path.exists(p):
            raise ValueError('Input folder from cellranger run (outs/)'
                             ' does not exist.')
    IOTools.touch_file(outfile)


def genClusterJobs():
    ''' Generate cluster jobs for each sample '''
    samples = pd.read_csv(PARAMS["input_samples"], sep='\t')
    infile = None
    samples.set_index("sample_id", inplace=True)

    for sample_name in samples.index:
        sample_dir = sample_name + ".sample.dir"
        out_sentinel = os.path.join(sample_dir, "sort.sentinel")
        yield(infile, out_sentinel)


# ########################################################################### #
# ############# Generate sorted bam files from cellranger bam ############### #
# ########################################################################### #

@follows(checkInputs)
@files(genClusterJobs)
def sortBam(infile, outfile):
    '''Sort bam file by cell barcodes'''

    sample_name = outfile.split("/")[0][:-len(".sample.dir")]

    # get path from input_samples.tsv
    samples = pd.read_csv(PARAMS["input_samples"], sep='\t')
    samples.set_index("sample_id", inplace=True)
    outfolder = samples.loc[sample_name, "path"]

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # cellranger outs directory contains input and output bam file for this task
    sort_outfile = os.path.join(outfolder, "cellsorted_possorted_genome_bam.bam")
    sort_infile = os.path.join(outfolder, "possorted_genome_bam.bam")
    log_file = outfile.replace("sentinel", "log")

    if not os.path.exists(sort_outfile):
        job_threads = PARAMS["sort_threads"]
        sort_threads = int(job_threads - 1)
        job_memory = PARAMS["sort_memory"]
        mem_ind = int(job_memory[:-len("M")])
        mem = int(mem_ind - mem_ind*0.1)
        sort_memory = str(mem) + "M"
        statement = '''samtools sort -t CB -O BAM
                       --threads %(sort_threads)s
                       -m %(sort_memory)s
                       -o %(sort_outfile)s %(sort_infile)s &> %(log_file)s
                    '''
    else:
        statement = '''echo 'Sorted file already present' > %(log_file)s
                    '''

    P.run(statement)

    # Create sentinel file
    IOTools.touch_file(outfile)



# ########################################################################### #
# ##################### Run velocyto on sorted bam ########################## #
# ########################################################################### #

@transform(sortBam,
           regex(r"(.*).sample.dir/sort.sentinel"),
           r"\1.sample.dir/run_velocity.sentinel")
def runVelocyto(infile, outfile):
    '''Run velocyto on barcode-sorted bam file. This task writes a loom file
       into the pipeline-run directory for each sample.'''

    sample_name = infile.split("/")[0][:-len(".sample.dir")]
    reference = os.path.join(str(PARAMS["velocyto_cellranger_anno"]),
                             "genes", "genes.gtf")
    outdir = os.path.dirname(outfile)

    samples = pd.read_csv(PARAMS["input_samples"], sep='\t')
    samples.set_index("sample_id", inplace=True)
    bcfile = samples.loc[sample_name, 'barcodes']
    bam_folder = samples.loc[sample_name, 'path']
    sort_outfile = os.path.join(bam_folder, "possorted_genome_bam.bam")
    log_file = outfile.replace("sentinel", "log")

    if "M" in PARAMS["velocyto_memory"]:
        job_memory = PARAMS["velocyto_memory"]

    job_threads = PARAMS["velocyto_threads"]
    statement = '''velocyto run --bcfile %(bcfile)s
                   --outputfolder %(outdir)s
                   --sampleid %(sample_name)s  -vvv
                   %(sort_outfile)s %(reference)s &> %(log_file)s
                '''
    P.run(statement)

    IOTools.touch_file(outfile)


# ---------------------------------------------------
# Generic pipeline tasks

@follows(runVelocyto)
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
