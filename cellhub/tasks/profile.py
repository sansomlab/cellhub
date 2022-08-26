'''
profile.py
==========

Parse the pipeline.log file from a cgat-core pipeline and
summarise information on task resource usage.

Usage
-----

After running a pipline, the resources used by the pipeline tasks can 
be summarised with the "cellhub pipeline_name profile" command, e.g.

.. code-block:: bash

    > cellhub cluster make full -v5 -p 100
    > cellhub cluster profile


Code
----

'''

import os
import re
import argparse
import pandas as pd
import logging
import sys
from shutil import which

import cgatcore.pipeline as P

# https://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script/34177358
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which
    return which(name) is not None

def setupLogger():

    L = logging.getLogger(__name__)
    log_handler = logging.StreamHandler(sys.stdout)
    log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
    log_handler.setLevel(logging.INFO)
    L.addHandler(log_handler)
    L.setLevel(logging.INFO)
    
    return L

def setupParser():

    parser = argparse.ArgumentParser()
    parser.add_argument("--log", default="pipeline.log", type=str,
                    help="File with reduced dimensions")
    parser.add_argument("--save-table", default=False, action="store_true",
                    help="Save the per-task tsv table.")
    
    
    return parser


def profile(log, save_table=False, show_fields=True):

    L = setupLogger()

    PARAMS = P.get_parameters()
    queue_manager = PARAMS["cluster_queue_manager"]
    
    if queue_manager not in ["slurm", "sge"]:
        raise ValueError("queue manager not supported")

    task_pattern = re.compile("[^{]*execution - {(.*)}$")

    tasks = {}
    n = 0
    with open (log, "r") as lf:

        for line in lf:

            if '"task"' in line and '"statement"' in line:

                t = task_pattern.match(line).group(1)
                t = t.replace("{","\{")
                t = t.replace("}","\}")
                t = t.replace(" null",'" null"')

                tasks[n] = eval("{" + t + "}")

                n+=1

    L.info("Parsing of log file complete")

    x = pd.DataFrame.from_dict(tasks, orient="index")

    if show_fields:
        L.info("Avaliable fields:")
        print(x.columns)

    if queue_manager == "slurm":
        x = x[["task",
            "NCPUS", "UserCPU","percent_cpu",
            "MaxVMSize","MaxRSS","MaxPages",
            "user_t","wall_t", "ExitCode"]]

        for mem_task in ["MaxVMSize", "MaxRSS"]:
            # we want this in GB
            x[mem_task] = x[mem_task] / 1E9 #1000000000000000

    elif queue_manager == "sge":

        x = x[["task",
            "slots", "percent_cpu",
            "max_vmem","max_rss","average_rss", "ru_nswap",
            "user_t","cpu_t", "wall_t", "exit_status"]]
        
        for mem_task in ["max_vmem", "max_rss","average_rss"]:
            # we want this in GB
            x[mem_task] = x[mem_task] / 1E6 #1000000000000000

    # I have no idea what the memory unit is for sge...!
    #
    #    for mem_task in ["max_vmem", "max_rss", "average_rss"]:
    #        x[mem_task] = x[mem_task] / 1000000000

    print("-----")
    L.info("Number of jobs per task:")
    print(x.groupby(['task']).size())

    print("-----")
    L.info("Average job resource usage:")
    print(x.groupby(['task']).mean().round(2))

    print("-----")
    L.info("Maximum job resource usage:")
    print(x.groupby(['task']).max().round(2))

    if save_table:
        x.to_csv("task.info.tsv.gz", sep="\t")

    L.info("Complete")


def main(argv=sys.argv):

    parser = setupParser()
    args = parser.parse_args()
        
    profile(args.log, save_table=args.save_table)


if __name__ == "__main__":
    sys.exit(main(sys.argv))