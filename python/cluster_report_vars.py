import os
import sys
import yaml
import logging
import argparse
from pathlib import Path
from types import SimpleNamespace
from cgatcore import pipeline as P
import cellhub.tasks.TASK as TASK

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help="infile")
parser.add_argument("--outfile",default=None, type=str,
                    help="outfile")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ############## Create outdir and set results file ######################### #
# ########################################################################### #

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline_cluster.yml" % os.path.splitext(__file__)[0],
     "../pipeline_cluster.yml",
     "pipeline_cluster.yml"])

# set the location of the code directory
PARAMS["cellhub_code_dir"] = Path(__file__).parents[1]

# get the task specification
spec, SPEC = TASK.get_vars(args.infile, args.outfile, PARAMS)
outfile_name = os.path.basename(args.outfile)

# initialise the namespace & alias the path function
x = SimpleNamespace()
p = os.path.join

# <---------------------------- base variables -----------------------------> #

x.outdir = spec.outdir
x.clusterDir = spec.cluster_dir
x.compDir = spec.component_dir
x.rdimsVisMethodShort = "umap"
x.clusterDirBaseName = os.path.basename(x.clusterDir)
x.nComponents = spec.components
x.resolution = spec.resolution
x.sample = Path(args.outfile).parts[0].split(".")[0]
x.sample = x.sample.replace("_", "\\_")

# <---------------------------- PARAMS variables -----------------------------> #

x.projectName = PARAMS["projectname"]
x.reportAuthor = PARAMS["author"]
x.cellhubDir = PARAMS["cellhub_code_dir"]
x.nnK = PARAMS["neighbors_n_neighbors"]
x.nnMethod = PARAMS["neighbors_method"]
x.nnMetric = PARAMS["neighbors_metric"]
x.threshUse = PARAMS["markers_min_fc"]
x.minPct = PARAMS["markers_min_pct"]

x.deTest = PARAMS["markers_test"]
x.clusteringAlgorithm= PARAMS["cluster_algorithm"]
x.reductionType = PARAMS["source_rdim_name"].replace("_", "\\_")
x.rdimsVisMethod = "umap.mindist_" + str(PARAMS["umap_mindist"])

# <------------------------------ path variables ---------------------------> #

x.singleRDir = p(x.compDir, "singleR.dir")
x.clusterGenesetsDir = p(x.clusterDir, "genesets.dir")
x.clusterMarkerDEPlotsDir = p(x.clusterDir, "marker.de.plots.dir")
x.clusterMarkerRdimsPlotsDir = p(x.clusterDir, "cluster.marker.rdims.plots.dir")
x.clusterMarkersDir = p(x.clusterDir, "markers.dir")
x.conditionGenesetsDir = p(x.clusterDir, "condition.genesets.dir")
x.conditionMarkerDEPlotsDir = p(x.clusterDir, "condition.marker.de.plots.dir")
x.conditionMarkersDir = p(x.clusterDir, "condition.markers.dir")
x.genelistsDir = p(x.clusterDir, "genelists.dir")
# x.knownmarkersDir = p(x.clusterDir, "known.markers.dir")
x.diffmapDir = p(x.clusterDir, "dm.visualisation.dir")
x.groupNumbersDir = p(x.clusterDir, "group.numbers.dir")
x.umapDir = p(x.compDir, "umap.dir")
x.rdimsVisClusterDir = p(x.clusterDir, "rdims.visualisation.dir")
x.rdimsVisFactorDir = p(x.compDir, "rdims.visualisation.dir")
x.rdimsVisSingleRDir = p(x.compDir, "singleR.dir", "rdims.visualisation.dir")
x.pagaDir = p(x.clusterDir,  "paga.dir")

# <------------------------------ blob variables ---------------------------> #

x.runName = x.nComponents + "\\_" + x.resolution
x.jobName = x.runName 

x.runDetails = ("no. components: " + str(x.nComponents) +
                ", cluster resolution: " + str(x.resolution) +
                ", cluster algorithm: " + str(x.clusteringAlgorithm) +
                ", de test: " + x.deTest)

# <-------------------------- conditional variables ------------------------> #

if PARAMS["markers_conserved"]:
    x.conservedFactor = PARAMS["markers_conserved_factor"]
    x.conservedFactor = x.conservedFactor.replace("_", "\\_")
else:
    x.conservedFactor = "None"

if PARAMS["markers_conserved_between"]:
    x.conservedBetweenFactor = PARAMS["markers_conserved_between_factor"]
    x.conservedBetweenFactor = x.conservedBetweenFactor.replace("_", "\\_")
else:
    x.conservedBetweenFactor = "None"


# <-------------------------- depreceated variables ------------------------> #

#x.sampleDir = spec.sample_dir
# x.phateDir = p(x.compDir, "phate.dir")
# x.velocityDir = p(x.compDir, "velocity.dir")
# x.qcMinGenes = PARAMS["qc_mingenes"]
# x.qcMaxMito = PARAMS["qc_maxpercentmito"]
# x.minCells = PARAMS["qc_mincells"]
# x.modelType = PARAMS["regress_modeluse"]
# x.cellCycle = PARAMS["regress_cellcycle"]
# x.sdCutOff = PARAMS["vargenes_sdcutoff"]
# x.latentVariables = PARAMS["regress_latentvars"].replace("_", "\\_")
# x.normalizationMethod = PARAMS["normalization_method"]
# x.nPositiveMarkers = PARAMS["exprsreport_n_positive"]
# x.nNegativeMarkers = PARAMS["exprsreport_n_negative"]


# <---------------------------- save the  variables ------------------------> #

with open(args.outfile, "w") as ofh:
    for command, value in x.__dict__.items():

        ofh.write("\\newcommand{\\" + command + "}{" + str(value) + "}\n")


L.info("Complete")
