import os
import sys


PROJECT_DIR = os.path.abspath(os.path.join(workflow.basedir, os.pardir))
sys.path.append(PROJECT_DIR)

from utils import parse2int, str2list

RSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/R/scripts"
PYSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/python"

# source
IN_ANNDATA = config["source"]["anndata"]
RDIM_NAME = config["source"]["rdim_name"]
HEATMAP_MAT = config["source"]["heatmap_matrix"]
CELLHUBAPI = config["source"]["cellhub"]
SINGLER_REFS = str2list(config["source"]["singler_refs"])  # NOTE: newly added

# runspecs
RDIMS_LST = [parse2int(x) for x in str2list(str(config["runspecs"]["n_components"]))]
MAX_RDIMS = max(RDIMS_LST)
RESOLUTION_LST = str2list(str(config["runspecs"]["cluster_resolutions"]))
PREDEFINED_CLUSTERS = config["runspecs"].get("predefined_clusters", None)

# run
GENE_IDS = "--gene_ids" if config["run"]["genesets"] else ""

# markers
CONSERVED = "--conserved" if config["markers"]["conserved"] else ""
CONSERVED_FACTOR = config["markers"]["conserved_factor"]
MARKERS_TEST = config["markers"]["test"]
MARKERS_PSEUDOCOUNT = config["markers"]["pseudocount"]


# NOTE: to be tested
def parse_subsetstat():
    import pandas as pd

    subset_factor = CONSERVED_FACTOR
    subset_stat = "--subset_factor=" + subset_factor
    levels_file = f"cluster.dir/metadata.dir/{CONSERVED_FACTOR}.levels"
    subset_levels = [x for x in pd.read_csv(levels_file, header=None)[0].values]
    return subset_levels


def get_clusters(ncomp, resolu):
    import pandas as pd

    cluster_ids = pd.read_csv(
        f"cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv",
        sep="\t",
        header=None,
    )[0].unique()
    return list(cluster_ids)


def all_marker_targets():
    from itertools import product

    targets = []
    for ncomp in RDIMS_LST:
        for resolu in RESOLUTION_LST:
            clusters = get_clusters(ncomp, resolu)
            for cluster, subset_level in product(clusters, SUBSET_LEVELS):
                targets.append(
                    f"cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/markers.dir/{cluster}.{subset_level}.markers.tsv.gz"
                )
    print(targets)
    return targets


SUBSET_LEVELS = parse_subsetstat() if CONSERVED else ["all"]


# neighbours
NEIGHBOUR_METHOD = config["neighbors"]["method"]
NEIGHBOUR_THREADS = config["neighbors"]["threads"]
NEIGHBOUR_K = config["neighbors"]["n_neighbors"]
NEIGHBOUR_METRIC = config["neighbors"]["metric"]
FULL_SPEED_MODE = "--fullspeed" if config["neighbors"]["full_speed"] else ""

# cluster
CLUSTER_ALGORITHM = config["cluster"]["algorithm"]

# UMAP
MIN_DIST = config["umap"]["mindist"]
MIN_DIST_LST = str2list(str(config["umap"]["mindists"]))

# plot
PDF = config["plot"]["pdf"]

GROUPS = str2list(config["plot"].get("groups", "cluster"))
if "cluster" not in GROUPS:
    GROUPS += ["cluster"]
SUBGROUPS = str2list(config["plot"].get("subgroup", None))
QCVARS = str2list(str(config["plot"]["qcvars"]))
RDIMCOLOURFACTORS = list(
    set([x for x in QCVARS + GROUPS + SUBGROUPS if x != "cluster"])
)

SHAPE = config["plot"].get("shape", None)
POINTALPHA = config["plot"]["pointalpha"]
POINTSIZE = config["plot"]["pointsize"]
POINTPCH = config["plot"]["pointpch"]

# summaries
SUMMARY_DICT = config["summaries"]


# Target files
TARGETS = (
    [
        "cluster.dir/preflight.log",  # preflight
        "cluster.dir/metadata.dir/metadata.tsv.gz",  # metadata
    ]
    + expand(
        "cluster.dir/loom.dir/{layer}.loom", layer=set([HEATMAP_MAT, "log1p"])
    )  # loom
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/{fl}",
        ncomp=RDIMS_LST,
        resolu=RESOLUTION_LST,
        fl=[
            "cluster_cell_counts.tsv",
            "cluster.dendrogram.png",
        ],  # scanpyCluster + clusterPostProcess + compareClusters
    )
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/clustree.{fmt}",
        ncomp=RDIMS_LST,
        fmt=["png", "pdf"],
    )  # clustree
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/{fl}",
        ncomp=RDIMS_LST,
        resolu=RESOLUTION_LST,
        fl=[
            "draw_graph_fa.paga.initialised.png",
            "draw_graph_fa.png",
            "paga_init_fa2.tsv.gz",
            "paga.png",
            "umap.paga.initialised.png",
            "umap.paga.init.tsv.gz",
        ],
    )  # PAGA
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/umap.dir/umap.{mindist}.tsv.gz",
        ncomp=RDIMS_LST,
        mindist=MIN_DIST_LST,
    )  # UMAP
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/rdims.visualisation.dir/UMAP.{fct}.png",
        ncomp=RDIMS_LST,
        fct=RDIMCOLOURFACTORS,
    )  # plotRdimsFactors
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/rdims.visualisation.dir/umap.mindist_{mindist}.cluster_id.png",
        ncomp=RDIMS_LST,
        resolu=RESOLUTION_LST,
        mindist=MIN_DIST_LST,
    )  # plotRdimsClusters
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/singleR.dir/UMAP.{ref}.pruned.labels.png",
        ncomp=RDIMS_LST,
        ref=SINGLER_REFS,
    )  # plotRdimsSingleR
    + expand(
        "cluster.dir/singleR.dir/{ref}.heatmap.png", ref=SINGLER_REFS
    )  # plotSingleR
    + ["cluster.dir/singleR.dir/summary.tex"]  # summariseSingleR
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir/number.plots.tex",
        ncomp=RDIMS_LST,
        resolu=RESOLUTION_LST,
    )  # summariseGroupNumbers
    + expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/stats.dir/{subset_level}.stats.tsv.gz",
        ncomp=RDIMS_LST,
        resolu=RESOLUTION_LST,
        subset_level=SUBSET_LEVELS,
    )  # clusterStats
    # + expand("{f}", f=all_marker_targets())
)


rule target:
    input:
        TARGETS,


rule preflight:
    input:
        IN_ANNDATA,
    output:
        "cluster.dir/preflight.log",
    log:
        "cluster.dir/preflight.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_preflight.py",
        rdim_name=RDIM_NAME,
        max_rdims=MAX_RDIMS,
        geneids=GENE_IDS,
        conserved=CONSERVED,
        conserved_factor=CONSERVED_FACTOR,
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --reduced_dims_name="{params.rdim_name}" \
            --max_reduced_dims="{params.max_rdims}" \
            {params.conserved} \
            {params.geneids} \
            --conserved_factor="{params.conserved_factor}" \
            &> "{log}"
        """


rule metadata:
    input:
        IN_ANNDATA,
    output:
        "cluster.dir/metadata.dir/metadata.tsv.gz",
    log:
        "cluster.dir/metadata.dir/metadata.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_metadata.py",
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
        IN_ANNDATA,
    output:
        "cluster.dir/loom.dir/{layer}.loom",
    log:
        "cluster.dir/loom.dir/{layer}.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_loom.py",
        layers="{layer}",
        outdir="cluster.dir/loom.dir",
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
        IN_ANNDATA,
    output:
        "cluster.dir/out.{ncomp}.comp.dir/neighbour.graph.h5ad",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/neighbour.graph.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_neighbor_graph.py",
        rdim_name=RDIM_NAME,
        ncomp="{ncomp}",
        method=NEIGHBOUR_METHOD,
        threads=NEIGHBOUR_THREADS,
        k=NEIGHBOUR_K,
        metric=NEIGHBOUR_METRIC,
        fullspeedmode=FULL_SPEED_MODE,
    threads: NEIGHBOUR_THREADS
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
        "cluster.dir/out.{ncomp}.comp.dir/neighbour.graph.h5ad",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/scanpy.clusters.tsv.gz",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/scanpy.clusters.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_cluster.py",
        ncomp="{ncomp}",
        algorithm=CLUSTER_ALGORITHM,
        resolution="{resolu}",
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir",
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --algorithm="{params.algorithm}" \
            --resolution="{params.resolution}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


# NOTE: PREDEFINED_CLUSTERS not sure
checkpoint clusterPostProcess:
    input:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/"
        + (PREDEFINED_CLUSTERS if PREDEFINED_CLUSTERS else "scanpy.clusters.tsv.gz"),
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_colors.tsv",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_cell_counts.tsv",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_postprocess.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_post_process.R",
        predefined=f"--predefined={PREDEFINED_CLUSTERS}" if PREDEFINED_CLUSTERS else "",
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/",
    shell:
        """
        Rscript "{params.script}" \
            --clusters="{input}" \
            {params.predefined} \
            --mincells=10 \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


def target_cluster_markers(wc):
    cb = checkpoints.clusterPostProcess.get(ncomp=wc.ncomp, resolu=wc.resolu)

    cluster_tsv = cb.output[0]

    if os.path.exists(cluster_tsv):
        ids = pd.read_csv(cluster_tsv, sep="\t", header=None, dtype=str)[0].tolist()
    else:
        raise FileNotFoundError(f"cluster ids not found at {cluster_tsv}.")

    return expand(
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/markers.dir/{cid}.{subset_level}.markers.tsv.gz",
        ncomp=wc.ncomp,
        resolu=wc.resolu,
        cid=cid,
    )


rule compareClusters:
    input:
        source_anndata=IN_ANNDATA,
        cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster.dendrogram.png",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/compare.clusters.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_compare.py",
        ncomp="{ncomp}",
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/",
        reductiontype=RDIM_NAME,
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input.source_anndata}" \
            --clusterids="{input.cluster_ids}" \
            --ncomp="{params.ncomp}" \
            --outdir="{params.outdir}" \
            --reduced_dims_name="{params.reductiontype}" \
            &> "{log}"
        """


rule clustTree:
    input:
        expand(
            "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
            ncomp=["{ncomp}"],
            resolu=RESOLUTION_LST,
        ),
    output:
        "cluster.dir/out.{ncomp}.comp.dir/clustree.png",
        "cluster.dir/out.{ncomp}.comp.dir/clustree.pdf",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/clustree.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_clustree.R",
        res_str=",".join(RESOLUTION_LST),
        id_files_str=lambda wildcards, input: ",".join(input),
        outdir="cluster.dir/out.{ncomp}.comp.dir/",
    shell:
        """
        Rscript "{params.script}" \
            --resolutions="{params.res_str}" \
            --clusteridfiles="{params.id_files_str}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule paga:
    input:
        neigh_anndata="cluster.dir/out.{ncomp}.comp.dir/neighbour.graph.h5ad",
        cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
        cluster_colours="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_colors.tsv",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/draw_graph_fa.paga.initialised.png",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/draw_graph_fa.png",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/paga_init_fa2.tsv.gz",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/paga.png",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/umap.paga.initialised.png",
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/umap.paga.init.tsv.gz",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/paga.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_paga.py",
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/paga.dir/",
    shell:
        """
        python "{params.script}" \
            --anndata="{input.neigh_anndata}" \
            --outdir="{params.outdir}" \
            --cluster_ids="{input.cluster_ids}" \
            --cluster_colors="{input.cluster_colours}" \
            &> "{log}"
        """


rule UMAP:
    input:
        "cluster.dir/out.{ncomp}.comp.dir/neighbour.graph.h5ad",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/umap.dir/umap.{mindist}.tsv.gz",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/umap.dir/umap.{mindist}.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_umap.py",
        mindist="{mindist}",
        outdir="cluster.dir/out.{ncomp}.comp.dir/umap.dir/",
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --mindist="{params.mindist}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


# NOTE: ggplot2 aes_string -> aes(!!sym)
# NOTE: possible to change the cluster_plot_rdims_factor as plotting per factor?
rule plotRdimsFactors:
    input:
        table=f"cluster.dir/out.{{ncomp}}.comp.dir/umap.dir/umap.{MIN_DIST}.tsv.gz",
        metadata="cluster.dir/metadata.dir/metadata.tsv.gz",
    output:
        expand(
            "cluster.dir/out.{{ncomp}}.comp.dir/rdims.visualisation.dir/UMAP.{fct}.png",
            fct=RDIMCOLOURFACTORS,
        )
        + [
            "cluster.dir/out.{ncomp}.comp.dir/rdims.visualisation.dir/UMAP.tex",
            "cluster.dir/out.{ncomp}.comp.dir/rdims.visualisation.dir/plot.rdims.factor.tex",
        ],
    log:
        "cluster.dir/out.{ncomp}.comp.dir/rdims.visualisation.dir/plot.rdims.factor.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_rdims_factor.R",
        colour_factor_arg="--colorfactors=" + ",".join(RDIMCOLOURFACTORS),
        shape_factor_arg=("--shapefactor=" + SHAPE) if SHAPE is not None else "",
        pointsize=POINTSIZE,
        pointalpha=POINTALPHA,
        pointpch=POINTPCH,
        pdf=PDF,
        outdir="cluster.dir/out.{ncomp}.comp.dir/rdims.visualisation.dir/",
    shell:
        """
        Rscript "{params.script}" \
            --table="{input.table}" \
            --metadata="{input.metadata}" \
            {params.colour_factor_arg} \
            {params.shape_factor_arg} \
            --pointsize="{params.pointsize}" \
            --pointalpha="{params.pointalpha}" \
            --pointpch="{params.pointpch}" \
            --pdf="{params.pdf}" \
            --outdir="{params.outdir}" \
            --plotdirvar=rdimsVisFactorDir \
            &> "{log}"
        """


rule plotRdimsClusters:
    input:
        table="cluster.dir/out.{ncomp}.comp.dir/umap.dir/umap.{mindist}.tsv.gz",
        cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/rdims.visualisation.dir/umap.mindist_{mindist}.cluster_id.png",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/rdims.visualisation.dir/plot.rdims.cluster.{mindist}.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_rdims_factor.R",
        umap_spec="umap.mindist_" + "{mindist}",
        shape_factor_arg=("--shapefactor=" + SHAPE) if SHAPE is not None else "",
        pointsize=POINTSIZE,
        pointalpha=POINTALPHA,
        pointpch=POINTPCH,
        pdf=PDF,
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/rdims.visualisation.dir/",
    shell:
        """
        Rscript "{params.script}" \
            --method="{params.umap_spec}" \
            --table="{input.table}" \
            --metadata="{input.cluster_ids}" \
            {params.shape_factor_arg} \
            --colorfactors=cluster_id \
            --pointsize="{params.pointsize}" \
            --pointalpha="{params.pointalpha}" \
            --pointpch="{params.pointpch}" \
            --pdf="{params.pdf}" \
            --outdir="{params.outdir}" \
            --plotdirvar=rdimsVisClusterDir \
            &> "{log}"
        """


# NOTE: added {ref} in log
rule plotRdimsSingleR:
    input:
        table=f"cluster.dir/out.{{ncomp}}.comp.dir/umap.dir/umap.{MIN_DIST}.tsv.gz",
        labels=os.path.join(CELLHUBAPI, "api", "singleR", "{ref}", "labels.tsv.gz"),
    output:
        "cluster.dir/out.{ncomp}.comp.dir/singleR.dir/UMAP.{ref}.pruned.labels.png",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/singleR.dir/rdims.plots.{ref}.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_rdims_factor.R",
        reference="{ref}",
        pointsize=POINTSIZE,
        pointalpha=POINTALPHA,
        pointpch=POINTPCH,
        pdf=PDF,
        outdir="cluster.dir/out.{ncomp}.comp.dir/singleR.dir",
    shell:
        """
        Rscript "{params.script}" \
            --table="{input.table}" \
            --metadata="{input.labels}" \
            --colorfactors=pruned.labels \
            --analysisname="{params.reference}" \
            --pointsize="{params.pointsize}" \
            --pointalpha="{params.pointalpha}" \
            --pointpch="{params.pointpch}" \
            --pdf="{params.pdf}" \
            --outdir="{params.outdir}" \
            --plotdirvar=rdimsVisClusterDir \
            &> "{log}"
        """


# NOTE: added {ref} in log
rule plotSingleR:
    input:
        metadata="cluster.dir/metadata.dir/metadata.tsv.gz",
        scores=os.path.join(CELLHUBAPI, "api", "singleR", "{ref}", "scores.tsv.gz"),
        labels=os.path.join(CELLHUBAPI, "api", "singleR", "{ref}", "labels.tsv.gz"),
    output:
        "cluster.dir/singleR.dir/{ref}.heatmap.png",
    log:
        "cluster.dir/singleR.dir/singleR.plots.{ref}.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_singleR_plots.R",
        reference="{ref}",
        outdir="cluster.dir/singleR.dir",
        pdf=PDF,
    shell:
        """
        Rscript "{params.script}" \
            --metadata="{input.metadata}" \
            --scores="{input.scores}" \
            --labels="{input.labels}" \
            --reference="{params.reference}" \
            --outdir="{params.outdir}" \
            --pdf="{params.pdf}" \
            &> "{log}"
        """


def summariseSingleR(
    singleR_path,
    ref_lst,
    out_path,
):
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


# NOTE: due to auto-formatter problem, we cannot put python script into `run` module.
rule summariseSingleR:
    input:
        expand("cluster.dir/singleR.dir/{ref}.heatmap.png", ref=SINGLER_REFS),
    output:
        "cluster.dir/singleR.dir/summary.tex",
    run:
        summariseSingleR(
            "cluster.dir/singleR.dir",
            SINGLER_REFS,
            "cluster.dir/singleR.dir/summary.tex",
        )


def populate_options(summary_key):
    options = []
    for k, v in SUMMARY_DICT[summary_key].items():
        if v == "None" or v == None or v == False or k == "title":
            pass
        elif v == True:
            options.append("--" + k)
        elif k in ["xlab", "ylab"]:
            options.append("--" + k + '="' + str(v) + '"')
        else:
            options.append("--" + k + '="' + str(v) + '"')
    return "\t".join(options)


# NOTE: cluster_plot_group_numbers.R aes_string -> aes(!!sym())
rule plotGroupNumbers:
    input:
        metadata="cluster.dir/metadata.dir/metadata.tsv.gz",
        cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir/{key}.data.tsv.gz",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir/plot.group.numbers.{key}.log",
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_group_numbers.R",
        key="{key}",
        options=lambda wc: populate_options(wc.key),
        outdir="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir",
    shell:
        """
        Rscript "{params.script}" \
            --metadata="{input.metadata}" \
            --clusters="{input.cluster_ids}" \
            --title="{params.key}" \
            {params.options} \
            --outdir="{params.outdir}" \
            --plotdirvar=groupNumbersDir \
            &> "{log}"
        """


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


# NOTE: split from the previous plotGroupNumbers
rule summariseGroupNumbers:
    input:
        expand(
            "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir/{key}.data.tsv.gz",
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            key=SUMMARY_DICT.keys(),
        ),
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/group.numbers.dir/number.plots.tex",
    params:
        outdir=lambda wc: f"cluster.dir/out.{wc.ncomp}.comp.dir/cluster.{wc.resolu}.dir/group.numbers.dir/",
    run:
        summariseGroupNumbers(SUMMARY_DICT, params.outdir)


rule clusterStats:
    input:
        anndata=IN_ANNDATA,
        cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv.gz",
    output:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/stats.dir/{subset_level}.stats.tsv.gz",
    log:
        "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/stats.dir/cluster.{subset_level}.stats.log",
    params:
        script=f"{PYSCRIPT_DIR}/cluster_stats.py",
        subset_stat="--subset_factor=" + CONSERVED_FACTOR if CONSERVED else "",
        subset_level="{subset_level}",
    shell:
        """
        python "{params.script}" \
            --anndata="{input.anndata}" \
            {params.subset_stat} \
            --subset_level="{params.subset_level}" \
            --clusterids="{input.cluster_ids}" \
            --outfile="{output}" \
            &> "{log}"
        """


# checkpoint checkClusters:
#     output:
#         cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv",
# # # NOTE: the 911 cluster - what is this?
# rule findMarkers:
#     input:
#         anndata=IN_ANNDATA,
#         cluster_ids="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/cluster_ids.tsv",
#     output:
#         "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/markers.dir/{cluster}.{subset_level}.markers.tsv.gz",
#     log:
#         "cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/markers.dir/{cluster}.{subset_level}.markers.log",
#     params:
#         script=f"{PYSCRIPT_DIR}/cluster_markers.py",
#         subset_stat="--subset_factor=" + CONSERVED_FACTOR if CONSERVED else "",
#         subset_level="{subset_level}",
#         cluster="{cluster}",
#         stats_file="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/stats.dir/{subset_level}.stats.tsv.gz",
#         sizes_file="cluster.dir/out.{ncomp}.comp.dir/cluster.{resolu}.dir/stats.dir/{subset_level}.sizes.tsv.gz",
#         markers_test=MARKERS_TEST,
#         markers_pseudocount=MARKERS_PSEUDOCOUNT,
#     shell:
#         """
#         if [ "{params.cluster}" != "911" ]; then
#             python "{params.script}" \
#                 --anndata="{input.anndata}" \
#                 {params.subset_stat} \
#                 --subset_level="{params.subset_level}" \
#                 --clusterids="{params.cluster}" \
#                 --group_means="{params.stats_file}" \
#                 --group_sizes="{params.sizes_file}" \
#                 --method="{params.markers_test}" \
#                 --pseudocount="{params.markers_pseudocount}" \
#                 --outfile="{output}" \
#                 &> {log}
#         fi
#         """
