import os
import sys
import yaml
import logging
import argparse
from pathlib import Path
import shutil
#from shutil import copyfile
from types import SimpleNamespace
from cgatcore import pipeline as P
#import cellhub.tasks.TASK as TASK

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
parser.add_argument("--latexvars", default=None, type=str,
                    help="the .sty file with the latex vars")
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
cellhub_code_dir = Path(__file__).parents[1]

source_dir = os.path.join(cellhub_code_dir,
                          "cellhub/reports/cluster_summary")

outfile_name = os.path.basename(args.outfile)

jobName = outfile_name[:-len(".pdf")]

outdir = os.path.dirname(args.outfile)
rundir = Path(outdir).parents[0]


latexVars = args.latexvars


# get the latex variables
s= ['''\\input %(latexVars)s''']
s.append('''\\def\\reportTitle{Cellhub cluster: summary report}''')
                
# get the intro
s.append('''\\input %(source_dir)s/introReport.tex''')
s.append('''\\input %(source_dir)s/introductionSection.tex''')
s.append('''\\input %(source_dir)s/taskSummary.tex''')

# add the section to visualise clusters and factors in reduced dimensions
# (plots made by tsne or umap)
s.append('''\\input %(source_dir)s/rdimsVisSection.tex''')

# singleR section
if(PARAMS["run_singleR"]):
    s.append('''\\input %(source_dir)s/singleRSection.tex''')

# add the section with plots of cell and gene numbers etc.
s.append('''\\input %(source_dir)s/numbersSection.tex''')

if(PARAMS["run_compare_clusters"]):
    s.append('''\\input %(source_dir)s/clusteringSection.tex''')

nresolutions = len(str(PARAMS["runspecs_cluster_resolutions"]).split(","))

if(nresolutions > 1):
    s.append('''\\input %(source_dir)s/clustree.tex''')

if(PARAMS["run_paga"]):
    s.append('''\\input %(source_dir)s/pagaSection.tex''')


# if(PARAMS["run_knownmarkers"]):
#    s.append('''\\input %(source_dir)s/knownmarkersSection.tex''')

s.append('''\\input %(source_dir)s/markerGenes.tex''')

if(PARAMS["run_top_marker_heatmap"]):
    s.append('''\\input %(source_dir)s/topMarkerHeatmap.tex''')

if(PARAMS["run_characterise_markers"]): # and not ...
    s.append('''\\input %(source_dir)s/markerGenesByCluster.tex''')

if(PARAMS["run_genesets"]):
    s.append('''\\input %(source_dir)s/genesetSection.tex''')

# When relevant, add section that compares
# two conditions within each cluster
if os.path.exists(
        os.path.join(rundir, "condition.markers.dir", "findMarkersBetweenConditions.sentinel")):

    wcc_section_name = "withinClusterComparisonSection.tex"
    s.append('''\\input %(source_dir)s/%(wcc_section_name)s''')

    if(PARAMS["run_genesets"]):
        s.append('''\\input %(source_dir)s/genesetBetweenSection.tex''')

s.append('''\\input %(cellhub_code_dir)s/cellhub/reports/latex/endmatter.tex''')

with open(args.outfile, "w") as out_file:

    out_file.write("\n".join(s) % locals())

L.info("Complete")
