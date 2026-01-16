import os
import textwrap
import pandas as pd

from types import SimpleNamespace
from pathlib import Path

from reports.template import template


def generate_report_vars(outfile, cluster, code_dir, ncomp, resolu):

    # ########################################################################### #
    # ############## Create outdir and set results file ######################### #
    # ########################################################################### #
    # load options from the config file
    PARAMS = cluster.params.copy()

    # set the location of the code directory
    PARAMS["cellhub_code_dir"] = os.path.abspath(code_dir)

    # initialise the namespace & alias the path function
    x = SimpleNamespace()

    # <---------------------------- base variables -----------------------------> #
    x.outdir = cluster.latex_dir(ncomp, resolu, abspath=False)
    x.clusterDir = cluster.cluster_dir(ncomp, resolu, abspath=False)
    x.compDir = cluster.comp_dir(ncomp, abspath=False)
    x.rdimsVisMethodShort = "umap"
    x.clusterDirBaseName = os.path.basename(x.clusterDir)

    x.nComponents = ncomp
    x.resolution = resolu
    # NOTE: to check - is x.sample used?
    x.sample = Path(x.outdir).parts[0].split(".")[0]
    x.sample = x.sample.replace("_", "\\_")

    # <---------------------------- PARAMS variables -----------------------------> #
    x.projectName = PARAMS["projectname"]
    x.reportAuthor = PARAMS["author"]
    x.cellhubDir = PARAMS["cellhub_code_dir"]

    x.nnK = PARAMS["neighbor_graph"]["n_neighbors"]
    x.nnMethod = PARAMS["neighbor_graph"]["method"]
    x.nnMetric = PARAMS["neighbor_graph"]["metric"]

    x.threshUse = PARAMS["markers"]["min_fc"]
    x.minPct = PARAMS["markers"]["min_pct"]
    x.deTest = PARAMS["markers"]["test"]

    x.clusteringAlgorithm = PARAMS["clustering"]["algorithm"]

    # NOTE: .replace should no longer be needed
    x.reductionType = PARAMS["dimension_reduction"]["rdim_name"]  # .replace("_", "\\_")

    x.rdimsVisMethod = "umap.mindist_" + str(PARAMS["plot"]["umap_mindist"])

    # <------------------------------ path variables ---------------------------> #
    x.umapDir = cluster.umap_dir(ncomp, abspath=False)
    x.rdimsVisFactorDir = cluster.rdims_vis_factor_dir(ncomp, abspath=False)
    x.groupNumbersDir = cluster.group_numbers_dir(ncomp, resolu, abspath=False)
    x.rdimsVisClusterDir = cluster.rdims_vis_cluster_dir(ncomp, resolu, abspath=False)
    x.clusterGenesetsDir = cluster.genesets_dir(ncomp, resolu, abspath=False)
    x.clusterMarkerDEPlotsDir = cluster.marker_de_plots_dir(
        ncomp, resolu, abspath=False
    )
    x.clusterMarkersDir = cluster.markers_dir(ncomp, resolu, abspath=False)
    x.pagaDir = cluster.paga_dir(ncomp, resolu, abspath=False)
    # x.singleRDir = p(x.compDir, "singleR.dir")
    # NOTE: may no longer be needed
    # x.clusterMarkerRdimsPlotsDir = cluster.marker_rdims_plots_dir(ncomp, resolu, abspath=False)
    # x.conditionGenesetsDir = p(x.clusterDir, "condition.genesets.dir")
    # x.conditionMarkerDEPlotsDir = p(x.clusterDir, "condition.marker.de.plots.dir")
    # x.conditionMarkersDir = p(x.clusterDir, "condition.markers.dir")
    # x.genelistsDir = p(x.clusterDir, "genelists.dir")
    # x.knownmarkersDir = p(x.clusterDir, "known.markers.dir")
    # x.diffmapDir = p(x.clusterDir, "dm.visualisation.dir")
    # x.rdimsVisSingleRDir = p(x.compDir, "singleR.dir", "rdims.visualisation.dir")

    # <------------------------------ blob variables ---------------------------> #
    x.runName = x.nComponents + "\\_" + x.resolution
    x.jobName = x.runName

    x.runDetails = (
        "no. components: "
        + str(x.nComponents)
        + ", cluster resolution: "
        + str(x.resolution)
        + ", cluster algorithm: "
        + str(x.clusteringAlgorithm)
        + ", de test: "
        + x.deTest
    )

    # <-------------------------- conditional variables ------------------------> #
    if PARAMS["markers"]["conserved"]:
        x.conservedFactor = PARAMS["markers"]["conserved_factor"]
        x.conservedFactor = x.conservedFactor.replace("_", "\\_")
    else:
        x.conservedFactor = "None"

    if PARAMS["markers"]["conserved_between"]:
        x.conservedBetweenFactor = PARAMS["markers"]["conserved_between_factor"]
        x.conservedBetweenFactor = x.conservedBetweenFactor.replace("_", "\\_")
    else:
        x.conservedBetweenFactor = "None"

    # <-------------------------- depreceated variables ------------------------> #
    # x.sampleDir = t.sample_dir
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
    with open(outfile, "w") as ofh:
        for command, value in x.__dict__.items():
            ofh.write("\\newcommand{\\" + command + "}{" + str(value) + "}\n")


def generate_summary_report(outfile, cluster, code_dir, ncomp, resolu):

    # ########################################################################### #
    # ############## Create outdir and set results file ######################### #
    # ########################################################################### #
    # load options from the config file
    PARAMS = cluster.params.copy()

    # set the location of the code directory
    source_dir = os.path.join(code_dir, "cellhub/reports/cluster_summary")

    latexvars = os.path.join(
        cluster.latex_dir(ncomp, resolu, abspath=False), "report.vars.sty"
    )

    outfile_name = os.path.basename(outfile)
    jobName = outfile_name[: -len(".tex")]

    outdir = os.path.dirname(outfile)
    rundir = os.path.abspath(os.path.join(outdir, os.pardir))

    # get the latex variables
    s = ["""\\input %(latexvars)s"""]
    s.append("""\\def\\reportTitle{Cellhub cluster: summary report}""")

    # get the intro
    s.append("""\\input %(source_dir)s/introReport.tex""")
    s.append("""\\input %(source_dir)s/introductionSection.tex""")
    s.append("""\\input %(source_dir)s/taskSummary.tex""")

    # add the section to visualise clusters and factors in reduced dimensions
    # (plots made by tsne or umap)
    s.append("""\\input %(source_dir)s/rdimsVisSection.tex""")

    # singleR section
    # if PARAMS["run_singleR"]:
    #     s.append("""\\input %(source_dir)s/singleRSection.tex""")

    # add the section with plots of cell and gene numbers etc.
    s.append("""\\input %(source_dir)s/numbersSection.tex""")

    if PARAMS["enable"]["compare_clusters"]:
        s.append("""\\input %(source_dir)s/clusteringSection.tex""")

    nresolutions = len(cluster.get_param("clustering", "resolutions", parse2list=True))
    if nresolutions > 1:
        s.append("""\\input %(source_dir)s/clustree.tex""")

    if PARAMS["enable"]["paga"]:
        s.append("""\\input %(source_dir)s/pagaSection.tex""")

    # if(PARAMS["run_knownmarkers"]):
    #    s.append('''\\input %(source_dir)s/knownmarkersSection.tex''')

    s.append("""\\input %(source_dir)s/markerGenes.tex""")

    if PARAMS["enable"]["top_marker_heatmap"]:
        s.append("""\\input %(source_dir)s/topMarkerHeatmap.tex""")

    if PARAMS["enable"]["characterise_markers"]:  # and not ...
        s.append("""\\input %(source_dir)s/markerGenesByCluster.tex""")

    if PARAMS["enable"]["genesets"]:
        s.append("""\\input %(source_dir)s/genesetSection.tex""")

    # When relevant, add section that compares
    # two conditions within each cluster
    if os.path.exists(
        os.path.join(
            rundir, "condition.markers.dir", "findMarkersBetweenConditions.sentinel"
        )
    ):
        wcc_section_name = "withinClusterComparisonSection.tex"
        s.append("""\\input %(source_dir)s/%(wcc_section_name)s""")
        if PARAMS["run_genesets"]:
            s.append("""\\input %(source_dir)s/genesetBetweenSection.tex""")

    s.append("""\\input %(code_dir)s/cellhub/reports/latex/endmatter.tex""")

    with open(outfile, "w") as out_file:
        out_file.write("\n".join(s) % locals() + "\n")


def _add_figure(plot_file=None, caption=None, width="1", height="0.9"):
    heatmap_fig = {"width": "1", "height": "0.9", "path": plot_file, "caption": caption}
    fig_tex = textwrap.dedent(template.figure % heatmap_fig)
    return fig_tex


def generate_marker_report(
    outfile, markers, latexvars, cluster, code_dir, ncomp, resolu
):

    # set the location of the code directory
    source_dir = os.path.join(code_dir, "cellhub/reports/cluster_marker")

    # not all clusters may have degenes
    markers = pd.read_csv(markers, sep="\t")
    clusters_with_markers = [x for x in markers.cluster.unique()]

    tex = []

    # <----------------------------- front matter ----------------------------> #
    tex.append("""\\input %(latexvars)s""")
    tex.append("""\\def\\reportTitle{CellHub cluster: marker report}""")
    tex.append("""\\input %(source_dir)s/clusterMarkerReport.tex""")

    # <----------------------------- overview plots ----------------------------> #
    tex.append(template.subsection % {"title": "overview plots"})

    tex.append(
        _add_figure(
            os.path.join(
                cluster.cluster_dir(ncomp, resolu),
                "rdims.visualisation.dir",
                "umap.mindist_"
                + str(cluster.get_param("plot", "umap_mindist"))
                + ".cluster_id",
            ),
            width="1",
            height=".9",
            caption="UMAP coloured by cluster  (resolution " + resolu + ")",
        )
    )

    tex.append("""\clearpage""")

    tmh = os.path.join(
        cluster.cluster_dir(ncomp, resolu), "markers.dir", "markers.summary.heatmap"
    )

    if os.path.exists(tmh + ".png"):
        tex.append(
            _add_figure(
                tmh,
                width="1",
                height="0.9",
                caption="Marker summary heatmap (resolution " + resolu + ")",
            )
        )

    tex.append("""\clearpage""")

    # <---------------------------- per cluster plots --------------------------> #
    for clust in clusters_with_markers:
        if str(clust) == "911":
            continue

        tex.append(
            template.subsection % {"title": "Markers for cluster: " + str(clust)}
        )

        fig = "heatmap"
        tex.append(
            _add_figure(
                os.path.join(
                    cluster.cluster_dir(ncomp, resolu),
                    "marker.plots.dir",
                    "cluster." + str(clust) + "." + fig,
                ),
                width="1",
                height="0.25",
                caption="cluster " + str(clust) + " " + fig,
            )
        )

        # fig = "dotplot"
        # tex.append(_add_figure(os.path.join(t.cluster_dir, "marker.plots.dir",
        #                             "cluster." + str(clust) + "." + fig),
        #             width = "1", height= "0.3",
        #             caption = "cluster " + str(clust) + " " +fig))

        fig = "rdims"
        tex.append(
            _add_figure(
                os.path.join(
                    cluster.cluster_dir(ncomp, resolu),
                    "marker.plots.dir",
                    "cluster." + str(clust) + "." + fig,
                ),
                width="1",
                height="0.5",
                caption="cluster " + str(clust) + " expression dot plot",
            )
        )

        # add the scatter plots.
        tex.append(
            _add_figure(
                os.path.join(
                    cluster.cluster_dir(ncomp, resolu),
                    "de.plots.dir",
                    "dePlots." + str(clust),
                ),
                width="1",
                height="0.9",
                caption="cluster " + str(clust) + " differential expression plots",
            )
        )

        # Add the violin plots
        fig = "violins"
        tex.append(
            _add_figure(
                os.path.join(
                    cluster.cluster_dir(ncomp, resolu),
                    "marker.plots.dir",
                    "cluster." + str(clust) + "." + fig,
                ),
                width="1",
                height="0.7",
                caption="cluster " + str(clust) + " " + fig,
            )
        )

    tex.append("""\\input %(code_dir)s/cellhub/reports/latex/endmatter.tex""")

    with open(outfile, "w") as out_file:
        out_file.write("\n".join(tex) % locals())
