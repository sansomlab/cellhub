'''
cluster.py
==========

Provides and extended setup class for pipeline_cluster.py

code
====

'''

import os
import pandas as pd
import cellhub.tasks.setup as setup

class setup(setup):

    def __init__(self, infile, outfile, PARAMS, 
                 memory="4G", cpu=1,
                 make_outdir=True):
    
        super().__init__(infile, outfile, PARAMS, 
                   memory=memory, cpu=cpu,
                   make_outdir=make_outdir,
                   expose_var=False)

        self.parts = outfile.split("/")
        self.nparts = len(self.parts)


        # get a list of the resolutions.
        if PARAMS["runspecs_cluster_resolutions"]:

            self.resolutions = [x.strip()
                                   for x in
                                    str(PARAMS["runspecs_cluster_resolutions"]).split(",")]

        else:
            self.resolutions = []

        if PARAMS["runspecs_predefined_clusters"]:
            self.resolutions.append("predefined")
            

        # if we are in the component directory
        if self.parts[0].endswith(".comp.dir"):

            self.component_dir = self.parts[0]
            self.components = self.parts[0].split(".")[1]
            
            self.neighbour_graph_path = os.path.join(self.component_dir,
                                            "neighbour.graph.h5ad")

            if os.path.exists(self.neighbour_graph_path):
               self.neighbour_graph_anndata = self.neighbour_graph_path

        # if we are in the cluster directory
        if self.nparts >= 2 and self.parts[1].startswith("cluster."):

            self.cluster_dir = os.path.join(self.parts[0], self.parts[1])
            self.resolution = self.parts[1][len("cluster."):-len(".dir")]
            self.cluster_ids = os.path.join(self.cluster_dir, "cluster_ids.tsv.gz")
            self.cluster_colors = os.path.join(self.cluster_dir, "cluster_colors.tsv")

            # get the list of clusters
            cluster_table = os.path.join(self.cluster_dir, "cluster_ids.tsv")
            
            if os.path.exists(cluster_table):
                clusters = pd.read_table(cluster_table, header=None)
                clusters = clusters[clusters.columns[0]].tolist()
                self.cluster_table = cluster_table
                self.clusters = [x for x in set(clusters)]
                self.nclusters = len(self.clusters)

        self.var = self.__dict__