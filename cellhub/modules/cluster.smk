import sys, os


PROJECT_DIR = os.path.abspath(os.path.join(workflow.basedir, os.pardir))
sys.path.append(PROJECT_DIR)


from utils import parse2int


RSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/R/scripts"
PYSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/python"
IN_ANNDATA = config["source"]["anndata"]

# source
RDIM_NAME = config["source"]["rdim_name"]
HEATMAP_MAT = config["source"]["heatmap_matrix"]

# runspecs
RDIMS_LST = [
    parse2int(x) for x in str(config["runspecs"]["n_components"]).strip().split(",")
]
MAX_RDIMS = max(RDIMS_LST)
RESOLUTION_LST = str(config["runspecs"]["cluster_resolutions"]).strip().split(",")
PREDEFINED_CLUSTERS = config["runspecs"].get("predefined_clusters", None)

# run
GENE_IDS = "--gene_ids" if config["run"]["genesets"] else ""

# markers
CONSERVED = "--conserved" if config["markers"]["conserved"] else ""
CONSERVED_FACTOR = config["markers"]["conserved_factor"]

# neighbours
NEIGHBOUR_METHOD = config["neighbors"]["method"]
NEIGHBOUR_THREADS = config["neighbors"]["threads"]
NEIGHBOUR_K = config["neighbors"]["n_neighbors"]
NEIGHBOUR_METRIC = config["neighbors"]["metric"]
FULL_SPEED_MODE = "--fullspeed" if config["neighbors"]["full_speed"] else ""

# cluster
CLUSTER_ALGORITHM = config["cluster"]["algorithm"]

# UMAP
MIN_DIST_LST = str(config["umap"]["mindists"]).strip().split(",")

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


# ??????????? PREDEFINED_CLUSTERS not sure
rule clusterPostProcess:
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
