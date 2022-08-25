import os

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


# def get(infile, outfile, PARAMS):
#     '''
#     Make a dictionary of commonly need paths and parameteres
#     '''

#     parts = outfile.split("/")
#     nparts = len(parts)

#     spec = { "sample_name": parts[0].split(".")[0],
#              "seurat_object": os.path.join(parts[0], "begin.rds"),
#              "outdir": os.path.dirname(outfile),
#              "indir": os.path.dirname(infile),
#              "resolutions" : [x.strip()
#                               for x in
#                               PARAMS["runspecs_cluster_resolutions"].split(",")]
#              }

#     # take care of making the output directory.
#     if not os.path.exists(spec["outdir"]):
#         os.mkdir(spec["outdir"])

#     if outfile.endswith(".sentinel"):
#         spec["log_file"] = outfile.replace(".sentinel", ".log")

#     # if we are in the component directory
#     if nparts >= 3:

#         spec["components"] = parts[1].split(".")[1]
#         graph_path = os.path.join(parts[0], parts[1], "graph.rds")

#         if os.path.exists(graph_path):
#             spec["graph"] =  graph_path

#     # if we are in the cluster directory
#     if nparts >= 4:

#         spec["resolution"] = parts[2][len("cluster."):-len(".dir")]

#     return(spec)