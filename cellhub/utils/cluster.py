import os
import textwrap
import pandas as pd

from types import SimpleNamespace
from pathlib import Path

from reports.template import template


def summariseSingleR(singleR_path, ref_lst, out_path):
    import os
    import textwrap
    from tasks.report import template as template

    singleR_umap_path = os.path.join(singleR_path, "umap")
    with open(out_path, "w") as tex:
        for reference in ref_lst:
            # heatmap
            tex.write(template.subsection % {"title": reference})
            tex.write("\n")
            heatmap_path = os.path.join(singleR_path, reference + ".heatmap")

            if os.path.exists(heatmap_path + ".png"):
                heatmap_fig = {
                    "width": "1",
                    "height": "0.9",
                    "path": heatmap_path,
                    "caption": "singleR predictions (" + reference + ")",
                }
                tex.write(textwrap.dedent(template.figure % heatmap_fig))
                tex.write("\n")

            umap_path = os.path.join(
                singleR_umap_path, "umap." + reference + ".pruned.labels"
            )

            if os.path.exists(umap_path + ".png"):
                umap_fig = {
                    "width": "1",
                    "height": "0.9",
                    "path": umap_path,
                    "caption": "pruned singleR predictions (" + reference + ")",
                }

                tex.write(textwrap.dedent(template.figure % umap_fig))
                tex.write("\n")


def summariseGroupNumbers(param_dict, outdir):
    import os
    import textwrap
    from tasks.report import template as template

    with open(os.path.join(outdir, "number.plots.tex"), "w") as tex:
        for fig in param_dict.keys():
            if "_" in param_dict[fig]["title"]:
                raise ValueError(
                    "Underscores are not allowed in the plot"
                    " titles (due to issues with latex..."
                )
            # Add the figures, one per subsection, escaping underscores.
            tex.write(template.subsection % {"title": param_dict[fig]["title"]})
            tex.write("\n")
            fig_path = os.path.join(outdir, fig)
            if os.path.exists(fig_path + ".png"):
                fig_spec = {
                    "width": "1",
                    "height": "0.9",
                    "path": fig_path,
                    "caption": param_dict[fig]["title"],
                }
                tex.write(textwrap.dedent(template.figure % fig_spec))
                tex.write("\n")


def generate_report_vars(outfile, path_dict, cluster, ncomp, resolu):

    # ########################################################################### #
    # ############## Create outdir and set results file ######################### #
    # ########################################################################### #
    # initialise the namespace & alias the path function
    x = SimpleNamespace()

    # <---------------------------- base variables -----------------------------> #
    print(path_dict)
    x.outdir = path_dict["latex"].format(ncomp=ncomp, resolu=resolu)
    x.clusterDir = path_dict["cluster"].format(ncomp=ncomp, resolu=resolu)
    x.compDir = path_dict["rdims"].format(ncomp=ncomp)
    x.rdimsVisMethodShort = "umap"
    # x.clusterDirBaseName = os.path.basename(x.clusterDir)

    x.nComponents = ncomp
    x.resolution = resolu
    # NOTE: to check - is x.sample used?
    x.sample = Path(x.outdir).parts[0].split(".")[0]
    x.sample = x.sample.replace("_", "\\_")

    # <---------------------------- PARAMS variables -----------------------------> #
    x.projectName = cluster.projectname
    x.reportAuthor = cluster.author
    x.cellhubDir = path_dict["cellhub_code_dir"]

    x.nnK = cluster.n_neigh
    x.nnMethod = cluster.neigh_method
    x.nnMetric = cluster.neigh_metric

    x.threshUse = cluster.min_fc
    x.minPct = cluster.min_pct
    x.deTest = cluster.test_method

    x.clusteringAlgorithm = cluster.clust_algo

    # NOTE: .replace should no longer be needed
    x.reductionType = cluster.rdim_name  # .replace("_", "\\_")

    x.rdimsVisMethod = "umap.mindist_" + str(cluster.main_mindist)

    # <------------------------------ path variables ---------------------------> #
    x.umapDir = path_dict["umap"].format(ncomp=ncomp)
    x.rdimsVisFactorDir = path_dict["rdims_factors"].format(ncomp=ncomp)
    x.groupNumbersDir = path_dict["group_numbers"].format(ncomp=ncomp, resolu=resolu)
    x.rdimsVisClusterDir = path_dict["rdims_clusters"].format(
        ncomp=ncomp, resolu=resolu
    )
    x.clusterGenesetsDir = path_dict["genesets"].format(ncomp=ncomp, resolu=resolu)
    x.clusterMarkerDEPlotsDir = path_dict["marker_de_plots"].format(
        ncomp=ncomp, resolu=resolu
    )
    x.clusterMarkersDir = path_dict["markers"].format(ncomp=ncomp, resolu=resolu)
    x.pagaDir = path_dict["paga"].format(ncomp=ncomp, resolu=resolu)
    x.singleRDir = path_dict["rdims_singler"].format(ncomp=ncomp)
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
    if cluster.conserved:
        x.conservedFactor = cluster.conserved_fact
        x.conservedFactor = x.conservedFactor.replace("_", "\\_")
    else:
        x.conservedFactor = "None"

    # if PARAMS["markers"]["conserved_between"]:
    #     x.conservedBetweenFactor = PARAMS["markers"]["conserved_between_factor"]
    #     x.conservedBetweenFactor = x.conservedBetweenFactor.replace("_", "\\_")
    # else:
    #     x.conservedBetweenFactor = "None"

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


def generate_summary_report(outfile, cluster, path_dir, ncomp, resolu):

    # ########################################################################### #
    # ############## Create outdir and set results file ######################### #
    # ########################################################################### #
    # set the location of the code directory
    source_dir = os.path.join(
        path_dir["cellhub_code_dir"], "cellhub/reports/cluster_summary"
    )

    latexvars = os.path.join(
        path_dir["latex"].format(ncomp=ncomp, resolu=resolu),
        "report.vars.sty",
    )

    # outfile_name = os.path.basename(outfile)
    # jobName = outfile_name[: -len(".tex")]

    outdir = os.path.dirname(outfile)
    rundir = os.path.abspath(os.path.join(outdir, os.pardir))

    # get the latex variables
    s = [f"\\input {latexvars}"]
    s.append("\\def\\reportTitle{Cellhub cluster: summary report}")

    # get the intro
    s.append(f"\\input {source_dir}/introReport.tex")
    s.append(f"\\input {source_dir}/introductionSection.tex")
    s.append(f"\\input {source_dir}/taskSummary.tex")

    # add the section to visualise clusters and factors in reduced dimensions
    # (plots made by tsne or umap)
    s.append(f"\\input {source_dir}/rdimsVisSection.tex")
    # singleR section
    if cluster.task_dict.get("singleR", False):
        s.append(f"\\input {source_dir}/singleRSection.tex")

    # add the section with plots of cell and gene numbers etc.
    s.append(f"\\input {source_dir}/numbersSection.tex")
    if cluster.task_dict.get("compare_clusters", False):
        s.append(f"\\input {source_dir}/clusteringSection.tex")

    nresolutions = len(cluster.clust_r_lst)
    if nresolutions > 1:
        s.append(f"\\input {source_dir}/clustree.tex")

    if cluster.task_dict.get("paga", False):
        s.append(f"\\input {source_dir}/pagaSection.tex")
    # if(PARAMS["run_knownmarkers"]):
    #    s.append('''\\input %(source_dir)s/knownmarkersSection.tex''')

    s.append(f"\\input {source_dir}/markerGenes.tex")

    if cluster.task_dict.get("top_marker_heatmap", False):
        s.append(f"\\input {source_dir}/topMarkerHeatmap.tex")

    if cluster.task_dict.get("characterise_markers", False):  # and not ...
        s.append(f"\\input {source_dir}/markerGenesByCluster.tex")

    if cluster.task_dict.get("genesets", False):
        s.append(f"\\input {source_dir}/genesetSection.tex")

    # When relevant, add section that compares
    # two conditions within each cluster
    # if os.path.exists(
    #     os.path.join(
    #         rundir, "condition.markers.dir", "findMarkersBetweenConditions.sentinel"
    #     )
    # ):
    #     wcc_section_name = "withinClusterComparisonSection.tex"
    #     s.append("""\\input %(source_dir)s/%(wcc_section_name)s""")
    #     if PARAMS["run_genesets"]:
    #         s.append("""\\input %(source_dir)s/genesetBetweenSection.tex""")

    s.append(
        f"\\input {path_dir['cellhub_code_dir']}/cellhub/reports/latex/endmatter.tex"
    )

    with open(outfile, "w") as out_file:
        out_file.write("\n".join(s) % locals() + "\n")


def _add_figure(plot_file=None, caption=None, width="1", height="0.9"):
    heatmap_fig = {
        "width": width,
        "height": height,
        "path": plot_file,
        "caption": caption,
    }
    fig_tex = textwrap.dedent(template.figure % heatmap_fig)
    return fig_tex


def generate_marker_report(
    outfile, markers, latexvars, path_dict, cluster, ncomp, resolu
):

    # set the location of the code directory
    source_dir = os.path.join(
        path_dict["cellhub_code_dir"], "cellhub/reports/cluster_marker"
    )

    # not all clusters may have degenes
    markers = pd.read_csv(markers, sep="\t")
    clusters_with_markers = [x for x in markers.cluster.unique()]

    tex = []

    # <----------------------------- front matter ----------------------------> #
    tex.append(f"\\input {latexvars}")
    tex.append("\\def\\reportTitle{CellHub cluster: marker report}")
    tex.append(f"\\input {source_dir}/clusterMarkerReport.tex")

    # <----------------------------- overview plots ----------------------------> #
    tex.append(template.subsection % {"title": "overview plots"})

    tex.append(
        _add_figure(
            os.path.join(
                path_dict["rdims_clusters"].format(ncomp=ncomp, resolu=resolu),
                "umap.mindist_" + str(cluster.main_mindist) + ".cluster_id",
            ),
            width="1",
            height="0.9",
            caption="UMAP coloured by cluster  (resolution " + resolu + ")",
        )
    )

    tex.append("""\clearpage""")

    tmh = os.path.join(
        path_dict["markers"].format(ncomp=ncomp, resolu=resolu),
        "markers.summary.heatmap",
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
                    path_dict["marker_plots"].format(ncomp=ncomp, resolu=resolu),
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
                    path_dict["marker_plots"].format(ncomp=ncomp, resolu=resolu),
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
                    path_dict["de_plots"].format(ncomp=ncomp, resolu=resolu),
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
                    path_dict["marker_plots"].format(ncomp=ncomp, resolu=resolu),
                    "cluster." + str(clust) + "." + fig,
                ),
                width="1",
                height="0.7",
                caption="cluster " + str(clust) + " " + fig,
            )
        )

    tex.append(
        f"\\input {path_dict['cellhub_code_dir']}/cellhub/reports/latex/endmatter.tex"
    )

    with open(outfile, "w") as out_file:
        out_file.write("\n".join(tex) % locals())
