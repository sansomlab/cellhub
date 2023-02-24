import os
import sys
import yaml
import logging
import argparse
import textwrap
from pathlib import Path
import shutil
import pandas as pd
#from shutil import copyfile
from types import SimpleNamespace
from cgatcore import pipeline as P
import cellhub.tasks.cluster as C
from cellhub.tasks.report import template as template

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
parser.add_argument("--markers", default=None, type=str,
                    help="the .tsv.gz file summarising the markers")
parser.add_argument("--outfile",default=None, type=str,
                    help="outfile")


args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ######################### function definitions ############################ #
# ########################################################################### #

def _add_figure(plot_file=None,
                caption=None,
                width="1",
                height="0.9"):

    heatmap_fig = {"width": "1", "height": "0.9",
                    "path": plot_file,
                    "caption": caption
    }

    fig_tex = textwrap.dedent(template.figure % heatmap_fig)

    return fig_tex
  


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
                          "cellhub/reports/cluster_marker")

t = C.setup(args.latexvars, args.outfile, PARAMS)

latexVars = args.latexvars

# not all clusters may have degenes
markers = pd.read_csv(args.markers, sep="\t")
markers.shape
#remove nas (this is added when looking for conserved markers)
markers = markers.loc[(~markers['p.adj'].isna()),:]
markers.shape
#filter significant genes
markers = markers.loc[markers['p.adj']<0.05]
clusters_with_markers = [x for x in markers.cluster.unique()]

L.info('cluster for report :')
L.info(clusters_with_markers)

tex = []

# <----------------------------- front matter ----------------------------> #

tex.append('''\\input %(latexVars)s''')
tex.append('''\\def\\reportTitle{CellHub cluster: marker report}''')
tex.append('''\\input %(source_dir)s/clusterMarkerReport.tex''')

# <----------------------------- overview plots ----------------------------> #

tex.append(template.subsection % {"title": "overview plots"})

tex.append(_add_figure(os.path.join(t.cluster_dir, "rdims.visualisation.dir",
                            "umap.mindist_" + str(PARAMS["umap_mindist"]) +\
                            ".cluster_id"),
                        width = "1", height= ".9",
            caption = "UMAP coloured by cluster  (resolution " + t.resolution + ")"))

tex.append('''\clearpage''')

tmh = os.path.join(t.cluster_dir, "markers.dir",
                    "markers.summary.heatmap")

if(os.path.exists(tmh + ".png")):
    tex.append(_add_figure(tmh,
                width = "1", height= "0.9",
                caption = "Marker summary heatmap (resolution " + t.resolution + ")"))

tex.append('''\clearpage''')

# <---------------------------- per cluster plots --------------------------> #

for clust in clusters_with_markers:

    if str(clust) == "911":
        continue

    tex.append(template.subsection % {"title": "Markers for cluster: " + str(clust)})

    fig = "heatmap"
    tex.append(_add_figure(os.path.join(t.cluster_dir, "marker.plots.dir",
                                "cluster." + str(clust) + "." + fig),
                width = "1", height= "0.25",
                caption = "cluster " + str(clust) + " " + fig))

    # fig = "dotplot"
    # tex.append(_add_figure(os.path.join(t.cluster_dir, "marker.plots.dir",
    #                             "cluster." + str(clust) + "." + fig),
    #             width = "1", height= "0.3",
    #             caption = "cluster " + str(clust) + " " +fig))

    fig = "rdims"
    tex.append(_add_figure(os.path.join(t.cluster_dir, "marker.plots.dir",
                                "cluster." + str(clust) + "." + fig),
                width = "1", height= "0.5",
                caption = "cluster " + str(clust) + " expression dot plot"))

    # add the scatter plots.
    tex.append(_add_figure(os.path.join(t.cluster_dir,
                                "de.plots.dir",
                                "dePlots." + str(clust)),
                width = "1", height= "0.9",
                caption = "cluster " + str(clust) + " differential expression plots"))

    # Add the violin plots
    fig = "violins"
    tex.append(_add_figure(os.path.join(t.cluster_dir, "marker.plots.dir",
                                "cluster." + str(clust) + "." + fig),
                width = "1", height= "0.7",
                caption = "cluster " + str(clust) + " " + fig))

tex.append('''\\input %(cellhub_code_dir)s/cellhub/reports/latex/endmatter.tex''')


with open(args.outfile, "w") as out_file:

    out_file.write("\n".join(tex) % locals())

L.info("Complete")
