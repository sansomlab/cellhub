import os
import math
import types
import pandas as pd

def get_resources(memory="4G", cpu=1, PARAMS={"resources_mempercore":False}):
    '''calculate the resource requirements and return a
       dictonary that can be used to update the local variables'''

    if not memory.endswith("G"):
        raise ValueError("Memory must be specified as XXG")

    gb_requested = int(memory[:-1])
    
    mpc = PARAMS["resources_mempercore"]

    mem_gb = int(math.ceil(gb_requested / float(cpu) ))

    if not mpc:
       
       ncpu = cpu
    
    else:
        # memory is requested via the core count
        # the number of GB per core is given in the resources_mempercore
        # (note that the schedule will simply ignore the job_memory)
        if not isinstance(mpc, int):
            raise ValueError("resources_mempercore must be False or integer")

        cpu_needed_for_mem_request = math.ceil(gb_requested / mpc)

        ncpu = max(cpu, cpu_needed_for_mem_request)


    spec = {"job_memory": str(mem_gb) + "G",
            "job_threads": ncpu,
            "r_memory": gb_requested * 1000}

    return(spec["job_threads"],
           spec["job_memory"],
           spec["r_memory"])


def get_vars(infile, outfile, PARAMS, make_outdir=True):
    '''
    Make a dictionary of commonly need paths and parameteres
    '''

    outfile = os.path.relpath(outfile)
    parts = outfile.split("/")

    nparts = len(parts)

    SPEC = { "outdir": os.path.dirname(outfile),
             "outname": os.path.basename(outfile),
             }

    # get a list of the resolutions.
    if PARAMS["runspecs_cluster_resolutions"]:

        SPEC["resolutions"] = [x.strip()
                               for x in
                               str(PARAMS["runspecs_cluster_resolutions"]).split(",")]

    else:
        SPEC["resolutions"] = []

    if PARAMS["runspecs_predefined_clusters"]:
        SPEC["resolutions"].append("predefined")


    if not infile is None:
        infile = os.path.relpath(infile)
        SPEC["indir"] = os.path.dirname(infile)
        SPEC["inname"] = os.path.basename(infile)

    if(make_outdir):

        # take care of making the output directory.
        if not os.path.exists(SPEC["outdir"]):
            os.makedirs(SPEC["outdir"])

    if outfile.endswith(".sentinel"):
        SPEC["log_file"] = outfile.replace(".sentinel", ".log")

    # if we are in the component directory
    if parts[0].endswith(".comp.dir"):

        SPEC["component_dir"] = parts[0]
        SPEC["components"] = parts[0].split(".")[1]

        neighbour_graph_path = os.path.join(SPEC["component_dir"],
                                    "neighbour.graph.h5ad")

        if os.path.exists(neighbour_graph_path):
            SPEC["neighbour_graph_anndata"] = neighbour_graph_path

    # if we are in the cluster directory
    if nparts >= 2 and parts[1].startswith("cluster."):

        SPEC["cluster_dir"] = os.path.join(parts[0], parts[1])
        SPEC["resolution"] = parts[1][len("cluster."):-len(".dir")]
        SPEC["cluster_ids"] = os.path.join(SPEC["cluster_dir"], "cluster_ids.tsv.gz")
        #SPEC["cluster_assignments"] = os.path.join(SPEC["cluster_dir"], "cluster_assignments.tsv.gz")
        SPEC["cluster_colors"] = os.path.join(SPEC["cluster_dir"], "cluster_colors.tsv")

        # get the list of clusters
        cluster_table = os.path.join(SPEC["cluster_dir"], "cluster_ids.tsv")
        if os.path.exists(cluster_table):
            clusters = pd.read_table(cluster_table, header=None)
            clusters = clusters[clusters.columns[0]].tolist()
            SPEC["cluster_table"] = cluster_table
            SPEC["clusters"] = [x for x in set(clusters)]
            SPEC["nclusters"] = len(SPEC["clusters"])

    return([types.SimpleNamespace(**SPEC), SPEC])
