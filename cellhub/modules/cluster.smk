import os
import sys
import pandas as pd

PROJECT_DIR = os.path.abspath(os.path.join(workflow.basedir, os.pardir))
sys.path.append(PROJECT_DIR)

from utils import parse2int, str2list, parse2float

RSCRIPT_DIR = os.path.join(PROJECT_DIR, "R", "scripts")
PYSCRIPT_DIR = os.path.join(PROJECT_DIR, "python")

ANNDATA_IN = config["anndata"]
OUT_DIR = config["out_dir"]

RDIMS_LST = [
    parse2int(ncomp, positive_only=True)
    for ncomp in str2list(str(config["dimension_reduction"]["n_components"]))
]
RESOLUTION_LST = [
    f"{parse2float(resolu, finite_only=True):.1f}"
    for resolu in str2list(config["clustering"]["resolutions"])
]
LAYER_LST = list(set(["log1p", config["plot"]["heatmap_matrix"]]))


def METADATA_DIR():
    return os.path.join(OUT_DIR, "metadata.dir")


def LOOM_DIR():
    return os.path.join(OUT_DIR, "loom.dir")


def RDIM_DIR(ncomp):
    return os.path.join(OUT_DIR, f"out.{ncomp}.comp.dir")


def UMAP_DIR(ncomp):
    return os.path.join(RDIM_DIR(ncomp), "umap.dir")


def CLUSTER_DIR(ncomp, resolu):
    return os.path.join(RDIM_DIR(ncomp), f"cluster.{resolu}.dir")


def STATS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "stats.dir")


def MARKERS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "markers.dir")


# NOTE: possible to combine MARKER PLOT FOLDERS?
def MARKER_PLOTS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "marker.plots.dir")


def DE_PLOTS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "de.plots.dir")


def MARKER_DE_PLOTS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "marker.de.plots.dir")


def GENESETS_DIR(ncomp, resolu):
    return os.path.join(CLUSTER_DIR(ncomp, resolu), "genesets.dir")


rule full:
    input:
        expand(
            os.path.join(RDIM_DIR("{ncomp}"), "clustree.png"),
            ncomp=RDIMS_LST,
        ),
        expand(
            os.path.join(
                MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.heatmap.png"
            ),
            ncomp=RDIMS_LST,
            resolu=RESOLUTION_LST,
        ),
        expand(
            os.path.join(
                DE_PLOTS_DIR("{ncomp}", "{resolu}"), "summarised_dePlots.sentinel"
            ),
            ncomp=RDIMS_LST,
            resolu=RESOLUTION_LST,
        ),
        expand(
            os.path.join(
                MARKER_PLOTS_DIR("{ncomp}", "{resolu}"),
                "summarised_markerPlots.sentinel",
            ),
            ncomp=RDIMS_LST,
            resolu=RESOLUTION_LST,
        ),
        expand(
            os.path.join(MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"), "deNumbers.png"),
            ncomp=RDIMS_LST,
            resolu=RESOLUTION_LST,
        ),
        expand(
            os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.xlsx"),
            ncomp=RDIMS_LST,
            resolu=RESOLUTION_LST,
        ),


# rule preflight:
#     input:
#         IN_ANNDATA,
#     output:
#         "cluster.dir/preflight.log",
#     log:
#         "cluster.dir/preflight.log",
#     params:
#         script=f"{PYSCRIPT_DIR}/cluster_preflight.py",
#         rdim_name=config["dimension_reduction"]["rdim_name"],
#         max_rdims=max(RDIMS_LST),
#         geneids=GENE_IDS,
#         conserved=CONSERVED,
#         conserved_factor=CONSERVED_FACTOR,
#     shell:
#         """
#         python "{params.script}" \
#             --anndata="{input}" \
#             --reduced_dims_name="{params.rdim_name}" \
#             --max_reduced_dims="{params.max_rdims}" \
#             {params.conserved} \
#             {params.geneids} \
#             --conserved_factor="{params.conserved_factor}" \
#             &> "{log}"
#         """


# NOTE: to implement conserved
CONSERVED = ""
CONSERVED_FACTOR = "stim"


rule metadata:
    input:
        ANNDATA_IN,
    output:
        os.path.join(METADATA_DIR(), "metadata.tsv.gz"),
    log:
        os.path.join(METADATA_DIR(), "metadata.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_metadata.py"),
        conserved=CONSERVED,
        conserved_factor=CONSERVED_FACTOR,
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input}" \
            {params.conserved} \
            --conserved_factor="{params.conserved_factor}" \
            --outfile="{output}" \
            &> "{log}"
        """


rule loom:
    input:
        ANNDATA_IN,
    output:
        os.path.join(LOOM_DIR(), "{layer}.loom"),
    log:
        os.path.join(LOOM_DIR(), "{layer}.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_loom.py"),
        layers="{layer}",
        outdir=LOOM_DIR(),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --layer="{params.layers}" \
            --loomdir="{params.outdir}" \
            &> "{log}"
        """


rule neighbourGraph:
    input:
        ANNDATA_IN,
    output:
        os.path.join(RDIM_DIR("{ncomp}"), "neighbour.graph.h5ad"),
    log:
        os.path.join(RDIM_DIR("{ncomp}"), "neighbour.graph.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_neighbor_graph.py"),
        rdim_name=config["dimension_reduction"]["rdim_name"],
        ncomp="{ncomp}",
        method=config["neighbor_graph"]["method"],
        threads=config["neighbor_graph"]["threads"],
        k=config["neighbor_graph"]["n_neighbors"],
        metric=config["neighbor_graph"]["metric"],
        fullspeedmode="--fullspeed" if config["neighbor_graph"]["full_speed"] else "",
    threads: config["neighbor_graph"]["threads"]
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input}" \
            --reduced_dims_name="{params.rdim_name}" \
            --outfile="{output}" \
            --ncomps="{params.ncomp}" \
            --method="{params.method}" \
            --threads="{params.threads}" \
            --k="{params.k}" \
            --metric="{params.metric}" \
            {params.fullspeedmode} \
            &> "{log}"
        """


rule scanpyCluster:
    input:
        os.path.join(RDIM_DIR("{ncomp}"), "neighbour.graph.h5ad"),
    output:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "scanpy.clusters.tsv.gz"),
    log:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "scanpy.clusters.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_cluster.py"),
        ncomp="{ncomp}",
        algorithm=config["clustering"]["algorithm"],
        resolution="{resolu}",
        outdir=CLUSTER_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --algorithm="{params.algorithm}" \
            --resolution="{params.resolution}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


# NOTE: to implement usage of pre-defined cluster column
checkpoint clusterPostProcess:
    """
    TODO: process predefined cluster
    """
    input:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "scanpy.clusters.tsv.gz"),
    output:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv"),
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_colors.tsv"),
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_cell_counts.tsv"),
    log:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_postprocess.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_post_process.R"),
        predefined="",  #f"--predefined={PREDEFINED_CLUSTERS}" if PREDEFINED_CLUSTERS else "",
        outdir=CLUSTER_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        Rscript "{params.script}" \
            --clusters="{input}" \
            {params.predefined} \
            --mincells=10 \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule compareClusters:
    input:
        anndata=ANNDATA_IN,
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
    output:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster.dendrogram.png"),
    log:
        os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "compare.clusters.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_compare.py"),
        ncomp="{ncomp}",
        outdir=CLUSTER_DIR("{ncomp}", "{resolu}"),
        reductiontype=config["dimension_reduction"]["rdim_name"],
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input.source_anndata}" \
            --clusterids="{input.cids}" \
            --ncomp="{params.ncomp}" \
            --outdir="{params.outdir}" \
            --reduced_dims_name="{params.reductiontype}" \
            &> "{log}"
        """


rule clustTree:
    input:
        expand(
            os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
            ncomp=["{ncomp}"],
            resolu=RESOLUTION_LST,
        ),
    output:
        os.path.join(RDIM_DIR("{ncomp}"), "clustree.png"),
        os.path.join(RDIM_DIR("{ncomp}"), "clustree.pdf"),
    log:
        os.path.join(RDIM_DIR("{ncomp}"), "clustree.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_clustree.R"),
        res_str=",".join(RESOLUTION_LST),
        id_files_str=lambda wc, input: ",".join(input),
        outdir=RDIM_DIR("{ncomp}"),
    shell:
        """
        Rscript "{params.script}" \
            --resolutions="{params.res_str}" \
            --clusteridfiles="{params.id_files_str}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule UMAP:
    input:
        os.path.join(RDIM_DIR("{ncomp}"), "neighbour.graph.h5ad"),
    output:
        os.path.join(UMAP_DIR("{ncomp}"), "umap.{mindist}.tsv.gz"),
    log:
        os.path.join(UMAP_DIR("{ncomp}"), "umap.{mindist}.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_umap.py"),
        mindist="{mindist}",
        outdir=UMAP_DIR("{ncomp}"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --mindist="{params.mindist}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


# NOTE: to implement subset_factor and subset_level
rule clusterStats:
    input:
        anndata=ANNDATA_IN,
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
    output:
        outfile1=os.path.join(
            STATS_DIR("{ncomp}", "{resolu}"), "{subset_level}.stats.tsv.gz"
        ),
        outfile2=os.path.join(
            STATS_DIR("{ncomp}", "{resolu}"), "{subset_level}.sizes.tsv.gz"
        ),
    log:
        os.path.join(
            STATS_DIR("{ncomp}", "{resolu}"), "cluster.{subset_level}.stats.log"
        ),
    params:
        script=f"{PYSCRIPT_DIR}/cluster_stats.py",
        subset_stat="",  # "--subset_factor=" + CONSERVED_FACTOR if CONSERVED else "",
        subset_level="{subset_level}",
    shell:
        """
        python "{params.script}" \
            --anndata="{input.anndata}" \
            {params.subset_stat} \
            --subset_level="{params.subset_level}" \
            --clusterids="{input.cids}" \
            --outfile="{output.outfile1}" \
            &> "{log}"
        """


def get_clusters(ncomp, resolu):
    cid_file = checkpoints.clusterPostProcess.get(ncomp=ncomp, resolu=resolu).output[0]
    with open(cid_file) as f:
        clusters = [
            line.strip() for line in f if line.strip() and (line.strip() != "911")
        ]
    return clusters


# NOTE: to implement subset_factor and subset_level
rule findMarkers:
    input:
        anndata=ANNDATA_IN,
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
        stats=os.path.join(
            STATS_DIR("{ncomp}", "{resolu}"), "{subset_level}.stats.tsv.gz"
        ),
        sizes=os.path.join(
            STATS_DIR("{ncomp}", "{resolu}"), "{subset_level}.sizes.tsv.gz"
        ),
    output:
        os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"),
            "{cluster}.{subset_level}.markers.tsv.gz",
        ),
    log:
        os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"),
            "{cluster}.{subset_level}.markers.log",
        ),
    params:
        script=f"{PYSCRIPT_DIR}/cluster_markers.py",
        subset_stat="",  # "--subset_factor=" + CONSERVED_FACTOR if CONSERVED else "",
        subset_level="{subset_level}",
        cluster="{cluster}",
        markers_test=config["markers"]["test"],
        markers_pseudocount=config["markers"]["pseudocount"],
    shell:
        """
        python "{params.script}" \
            --anndata="{input.anndata}" \
            {params.subset_stat} \
            --subset_level="{params.subset_level}" \
            --clusterids="{input.cids}" \
            --cluster="{params.cluster}" \
            --group_means="{input.stats}" \
            --group_sizes="{input.sizes}" \
            --method="{params.markers_test}" \
            --pseudocount="{params.markers_pseudocount}" \
            --outfile="{output}" \
            &> {log}
        """


def get_cluster_marker_files(ncomp, resolu):
    return [
        os.path.join(
            MARKERS_DIR(ncomp, resolu),
            f"{cluster}.{subset_level}.markers.tsv.gz",
        )
        for cluster in get_clusters(ncomp, resolu)
        for subset_level in ["all"]
    ]


checkpoint summariseMarkers:
    input:
        metadata=os.path.join(OUT_DIR, "metadata.dir", "metadata.tsv.gz"),
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
        cluster_markers=lambda wc: get_cluster_marker_files(wc.ncomp, wc.resolu),
    output:
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"),
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.xlsx"),
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.stats.tsv"),
    log:
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers_summary.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_summarise_markers.R"),
        markers_str=lambda wc: ",".join(get_cluster_marker_files(wc.ncomp, wc.resolu)),
        min_pct=config["markers"]["min_pct"],
        min_fc=config["markers"]["min_fc"],
        outdir=MARKERS_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        Rscript "{params.script}" \
            --marker_files="{params.markers_str}" \
            --minpct="{params.min_pct}" \
            --minfc="{params.min_fc}" \
            --clusterids="{input.cids}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule topMarkerHeatmap:
    input:
        loom_file=os.path.join(
            OUT_DIR, "loom.dir", f"{config['plot']['heatmap_matrix']}.loom"
        ),
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
        metadata=os.path.join(OUT_DIR, "metadata.dir", "metadata.tsv.gz"),
        marker_table=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
        ),
    output:
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.heatmap.png"),
    log:
        os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "topMarkerHeatmap.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_top_marker_heatmap.R"),
        matrix_loc="matrix",
        scale=(config["plot"]["heatmap_matrix"] == "log1p"),
        pdf=config["plot"]["pdf"],
        subgroup=(
            f"--subgroup={config['plot']['subgroup']}"
            if config["plot"]["subgroup"]
            else ""
        ),
        outdir=MARKERS_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        Rscript "{params.script}" \
            --loom="{input.loom_file}" \
            --clusterids="{input.cids}" \
            --metadata="{input.metadata}" \
            --matrix_loc="{params.matrix_loc}" \
            --scale="{params.scale}" \
            --markers="{input.marker_table}" \
            --pdf="{params.pdf}" \
            {params.subgroup} \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


# NOTE: parallelised
rule dePlots:
    input:
        marker_table=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
        ),
    output:
        os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.png"),
        os.path.join(
            DE_PLOTS_DIR("{ncomp}", "{resolu}"), "characterise.degenes.{cluster}.tex"
        ),
    log:
        os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_de_plots.R"),
        cluster="{cluster}",
        outdir=DE_PLOTS_DIR("{ncomp}", "{resolu}"),
        pdf=config["plot"]["pdf"],
    shell:
        """
        Rscript "{params.script}" \
            --degenes="{input.marker_table}" \
            --cluster="{params.cluster}" \
            --outdir="{params.outdir}" \
            --pdf="{params.pdf}" \
            --plotdirvar=clusterMarkerDEPlotsDir \
            &> "{log}"
        """


rule summarise_dePlots:
    input:
        expand(
            os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.png"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=lambda wc: get_clusters(wc.ncomp, wc.resolu),
        ),
        expand(
            os.path.join(
                DE_PLOTS_DIR("{ncomp}", "{resolu}"),
                "characterise.degenes.{cluster}.tex",
            ),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=lambda wc: get_clusters(wc.ncomp, wc.resolu),
        ),
    output:
        os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "summarised_dePlots.sentinel"),
    shell:
        """
        touch "{output}"
        """


def get_clusters_with_marker(ncomp, resolu):
    markerfile = checkpoints.summariseMarkers.get(ncomp=ncomp, resolu=resolu).output[0]
    markers = pd.read_csv(markerfile, sep="\t")
    clusters_with_markers = list(markers["cluster"].unique())
    assert (
        "911" not in clusters_with_markers
    ), "cluster 911 should not used for markerPlots."

    return clusters_with_markers


# NOTE: parallelised
# NOTE: markers in CGAT version line 1244?
# NOTE: group_opt not implemented
# NOTE: violinplot not visible when few clusters
rule markerPlots:
    input:
        marker_table=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
        ),
        loom_file=os.path.join(LOOM_DIR(), "log1p.loom"),
        metadata_file=os.path.join(METADATA_DIR(), "metadata.tsv.gz"),
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
        rdims_table=os.path.join(
            UMAP_DIR("{ncomp}"), f"umap.{config['plot']['umap_mindist']}.tsv.gz"
        ),
    output:
        rdims=os.path.join(
            MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.rdims.png"
        ),
        violins=os.path.join(
            MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.violins.png"
        ),
        heatmap=os.path.join(
            MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.heatmap.png"
        ),
    log:
        os.path.join(
            MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "marker.plots.{cluster}.log"
        ),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_marker_plots.R"),
        scaled_loom=(
            f"--scaled_loom={os.path.join(LOOM_DIR(), 'X.loom')}"
            if config["plot"]["heatmap_matrix"] == "X"
            else ""
        ),
        rdims_table=os.path.join(
            UMAP_DIR("{ncomp}"), f"umap.{config['plot']['umap_mindist']}.tsv.gz"
        ),
        cluster="{cluster}",
        outdir=MARKER_PLOTS_DIR("{ncomp}", "{resolu}"),
        group_opt="",
        pdf=config["plot"]["pdf"],
    shell:
        """
        Rscript "{params.script}" \
            --markers="{input.marker_table}" \
            --loom="{input.loom_file}" \
            {params.scaled_loom} \
            --metadata="{input.metadata_file}" \
            --clusterids="{input.cids}" \
            --rdimstable="{input.rdims_table}" \
            --cluster="{params.cluster}" \
            --outdir="{params.outdir}" \
            {params.group_opt} \
            --pdf="{params.pdf}" \
            &> "{log}"
        """


rule summarise_markerPlots:
    input:
        expand(
            os.path.join(
                MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.{plot}.png"
            ),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
            plot=["rdims", "violins", "heatmap"],
        ),
    output:
        os.path.join(
            MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "summarised_markerPlots.sentinel"
        ),
    shell:
        """
        touch "{output}"
        """


# NOTE: minfc and minpadj hard coded?
rule plotMarkerNumbers:
    input:
        marker_table=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
        ),
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
    output:
        os.path.join(MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"), "deNumbers.png"),
    log:
        os.path.join(
            MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"), "plotMarkerNumbers.log"
        ),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_marker_numbers.R"),
        outdir=MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        Rscript "{params.script}" \
            --degenes="{input.marker_table}" \
            --clusterids="{input.cids}" \
            --outdir="{params.outdir}" \
            --minfc=2 \
            --minpadj=0.05 \
            --plotdirvar=clusterMarkerDEPlotsDir \
            &> "{log}"
        """


# NOTE: re-written
def parseGMTfiles(config_contents):
    """
    Helper function for parsing the lists of GMT files
    """
    all_files = []
    for gmt_dict in config_contents:
        if gmt_dict is not None:
            for gmt_file in gmt_dict.values():
                if gmt_file is not None:
                    all_files += [gmt_file]
    if len(all_files) == 0:
        all_files = "none"
    else:
        all_files = ",".join(all_files)
    return all_files


# NOTE: re-written
def parseGMTnames(config_contents):
    """
    Helper function for parsing the lists of GMT files
    """
    all_names = []
    for gmt_dict in config_contents:
        if gmt_dict is not None:
            for gmt_name in gmt_dict.keys():
                if gmt_name is not None:
                    all_names += [gmt_name]
    if len(all_names) == 0:
        all_names = "none"
    else:
        all_names = ",".join(all_names)
    return all_names


# NOTE: replaced the dependence of CellHub API by actual annotation files
rule genesetAnalysis:
    input:
        marker_table=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
        ),
    output:
        expand(
            os.path.join(
                GENESETS_DIR("{ncomp}", "{resolu}"),
                "genesets.{cluster}.{files}.tsv.gz",
            ),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=["{cluster}"],
            files=[
                "GO.BP",
                "GO.CC",
                "GO.MF",
                "KEGG",
                "msigdb_biocarta",
                "msigdb_reactome",
            ],
        ),
    log:
        os.path.join(
            GENESETS_DIR("{ncomp}", "{resolu}"), "geneset.analysis.{cluster}.log"
        ),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_geneset_analysis.R"),
        cluster="{cluster}",
        universe=os.path.join(
            MARKERS_DIR("{ncomp}", "{resolu}"), "{cluster}.universe.tsv.gz"
        ),
        species=config["geneset"]["species"],
        ensembl=config["geneset"]["cellhub_ensembl_annotations"],
        kegg=config["geneset"]["cellhub_kegg_pathways"],
        gmt_names=parseGMTnames(
            [config["gmt_celltype_files"], config["gmt_pathway_files"]]
        ),
        gmt_files=parseGMTfiles(
            [config["gmt_celltype_files"], config["gmt_pathway_files"]]
        ),
        adjpthreshold=config["geneset"]["marker_adjpthreshold"],
        outdir=GENESETS_DIR("{ncomp}", "{resolu}"),
    shell:
        """
        Rscript "{params.script}" \
            --markers="{input.marker_table}" \
            --universe="{params.universe}" \
            --species="{params.species}" \
            --annotation="{params.ensembl}" \
            --kegg_pathways="{params.kegg}" \
            --gmt_names="{params.gmt_names}" \
            --gmt_files="{params.gmt_files}" \
            --cluster="{params.cluster}" \
            --adjpthreshold="{params.adjpthreshold}" \
            --direction=positive \
            --outdir="{params.outdir}" \
            &> "{log}"
        for output in "{output}"; do
            touch $output
        done
        """


# NOTE: possible to have clusters with no gene sets?
rule summariseGenesetAnalysis:
    input:
        lambda wc: [
            os.path.join(
                GENESETS_DIR(wc.ncomp, wc.resolu), f"genesets.{cluster}.{file}.tsv.gz"
            )
            for file in [
                "GO.BP",
                "GO.CC",
                "GO.MF",
                "KEGG",
                "msigdb_biocarta",
                "msigdb_reactome",
            ]
            for cluster in get_clusters_with_marker(wc.ncomp, wc.resolu)
        ],
        cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
    output:
        os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.xlsx"),
        os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.table.tex"),
        os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.figure.tex"),
    log:
        os.path.join(
            GENESETS_DIR("{ncomp}", "{resolu}"), "summarise.geneset.analysis.log"
        ),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_geneset_summary.R"),
        geneset_dir=GENESETS_DIR("{ncomp}", "{resolu}"),
        gmt_names=parseGMTnames(
            [config["gmt_celltype_files"], config["gmt_pathway_files"]]
        ),
        show_detailed=config["geneset"]["show_detailed"],
        min_genes=config["geneset"]["min_fg_genes"],
        pvalue_threshold=config["geneset"]["pvalue_threshold"],
        padjust_method=config["geneset"]["padjust_method"],
        use_adjusted=config["geneset"]["use_adjusted_pvalues"],
        min_odds_ratio=config["geneset"]["min_odds_ratio"],
        show_common=config["geneset"]["show_common"],
        out_prefix=os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets"),
        pdf=config["plot"]["pdf"],
    shell:
        """
        Rscript "{params.script}" \
            --genesetdir="{params.geneset_dir}" \
            --gmt_names="{params.gmt_names}" \
            --show_detailed="{params.show_detailed}" \
            --clusters="{input.cids}" \
            --mingenes="{params.min_genes}" \
            --pvaluethreshold="{params.pvalue_threshold}" \
            --padjustmethod="{params.padjust_method}" \
            --useadjusted="{params.use_adjusted}" \
            --minoddsratio="{params.min_odds_ratio}" \
            --showcommon="{params.show_common}" \
            --outprefix="{params.out_prefix}" \
            --prefix="genesets" \
            --plotdirvar="clusterGenesetsDir" \
            --pdf="{params.pdf}" \
            &> "{log}"
        """
