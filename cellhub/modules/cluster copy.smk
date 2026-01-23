import os
import sys
import pandas as pd

CELLHUB_CODE_DIR = os.path.abspath(os.path.join(workflow.basedir, os.pardir))
sys.path.append(CELLHUB_CODE_DIR)

from parser.parse_args import parse2int, str2list, parse2float
from parser.cluster_setup import ClusterSetup

cluster = ClusterSetup(config)

RSCRIPT_DIR = os.path.join(CELLHUB_CODE_DIR, "R", "scripts")
PYSCRIPT_DIR = os.path.join(CELLHUB_CODE_DIR, "python")
LATEX_DIR = os.path.join(CELLHUB_CODE_DIR, "cellhub", "reports")


# <------------------------- Utility Functions -------------------------> #
# NOTE: to be called in Snakemake after metadata task is done
def get_conserved_levels(wc):
    levels_path = cluster.path_tpl["metadata"]["outputs"].get("levels")
    if levels_path is None:
        return ["all"]
    ck = checkpoints.metadata.get(**wc)
    return pd.read_csv(levels_path, header=None)[0].tolist()


def get_valid_clusters(ncomp, resolu):
    ck = checkpoints.cluster_post_process.get(ncomp=ncomp, resolu=resolu)
    cid_file = ck.output["cids_uq"]
    return [x for x in pd.read_csv(cid_file, header=None)[0] if x != 911]


def get_cluster_marker_files(ncomp, resolu):

    return [
        cluster.path_tpl["find_markers"]["output"].format(
            ncomp="{ncomp}", resolu="{resolu}", level="{level}", cluster=cluster
        )
        for cluster in get_valid_clusters(ncomp, resolu)
    ]


# <------------------------- Exit Rules -------------------------> #
rule core:
    input:
        cluster.get_task_summary_outfile(),
        cluster.get_preflight_log(),
        cluster.get_metadata_outfiles(),
        cluster.get_clustree_outfiles(),
        cluster.get_rdims_vis_factors_outfiles(),
        cluster.get_rdims_vis_clusters_outfiles(),
        cluster.get_summarise_group_numbers_outfiles(),
        lambda wc: cluster.get_cluster_stats_outfiles(get_conserved_levels(wc)),


rule optional:
    input:
        cluster.get_compare_clusters_outfiles(),
        cluster.get_paga_outfiles(),
        cluster.get_summarise_singler_outfiles(),
        lambda wc: cluster.get_find_markers_outfiles(
            get_valid_clusters(wc), get_conserved_levels(wc)
        ),


rule full:
    input:
        rules.core.input,
        rules.optional.input,


# <------------------------- Process Rules -------------------------> #
rule task_summary:
    output:
        cluster.path_tpl["task_summary"]["output"],
    run:
        tasks = cluster.task_dict.keys()
        runs = cluster.task_dict.values()
        tab = pd.DataFrame.from_dict({"task": tasks, "run": runs})
        tab.to_latex(buf=output[0], index=False, escape=True)


rule preflight:
    input:
        cluster.anndata,
    output:
        cluster.path_tpl["preflight"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_preflight.py"),
        rdim_name=cluster.rdim_name,
        max_rdims=cluster.max_rdims,
        conserved=cluster.conserved_arg,
        geneids=cluster.preflight_geneids,
        conserved_factor=cluster.conserved_fact,
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --reduced_dims_name="{params.rdim_name}" \
            --max_reduced_dims="{params.max_rdims}" \
            {params.conserved} \
            {params.geneids} \
            --conserved_factor="{params.conserved_factor}" \
            &> "{output[0]}"
        """


checkpoint metadata:
    input:
        cluster.anndata,
    output:
        list(cluster.path_tpl["metadata"]["outputs"].values()),
    log:
        cluster.path_tpl["metadata"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_metadata.py"),
        conserved=cluster.conserved_arg,
        conserved_factor=cluster.conserved_fact,
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input}" \
            {params.conserved} \
            --conserved_factor="{params.conserved_factor}" \
            --outfile="{output[0]}" \
            &> "{log}"
        """


rule loom:
    input:
        cluster.anndata,
    output:
        cluster.path_tpl["loom"]["output"],
    log:
        cluster.path_tpl["loom"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_loom.py"),
        layers="{layer}",
        outdir=cluster.path_tpl["loom"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --layer="{params.layers}" \
            --loomdir="{params.outdir}" \
            &> "{log}"
        """


rule neighbour_graph:
    input:
        cluster.anndata,
    output:
        cluster.path_tpl["neighbour_graph"]["output"],
    log:
        cluster.path_tpl["neighbour_graph"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_neighbor_graph.py"),
        rdim_name=cluster.rdim_name,
        ncomp="{ncomp}",
        method=cluster.neigh_method,
        k=cluster.n_neigh,
        metric=cluster.neigh_metric,
        columns_to_keep=(
            f"--keep_obs {cluster.predef_clust_col}"
            if cluster.predef_clust_col
            else ""
        ),
        threads=cluster.hnsw_threads,
        fullspeedmode="--fullspeed" if cluster.hnsw_fullspeed else "",
    threads: cluster.hnsw_threads
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input}" \
            --reduced_dims_name="{params.rdim_name}" \
            --outfile="{output}" \
            --ncomps="{params.ncomp}" \
            --method="{params.method}" \
            --k="{params.k}" \
            --metric="{params.metric}" \
            {params.columns_to_keep} \
            --threads="{params.threads}" \
            {params.fullspeedmode} \
            &> "{log}"
        """


rule scanpy_cluster:
    input:
        rules.neighbour_graph.output,
    output:
        cluster.path_tpl["scanpy_cluster"]["output"],
    log:
        cluster.path_tpl["scanpy_cluster"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_cluster.py"),
        ncomp="{ncomp}",
        resolution="{resolu}",
        algorithm=cluster.clust_algo,
        cluster_colname=(
            f"--cluster_colname {cluster.predef_clust_col}"
            if cluster.clust_algo == "predefined"
            else ""
        ),
        outdir=cluster.path_tpl["scanpy_cluster"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --algorithm="{params.algorithm}" \
            --resolution="{params.resolution}" \
            {params.cluster_colname} \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


checkpoint cluster_post_process:
    input:
        rules.scanpy_cluster.output,
    output:
        cids_uq = cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_uq"],
        cids_full = cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
        ccolors = cluster.path_tpl["cluster_postprocess"]["outputs"]["ccolors"],
        
    log:
        cluster.path_tpl["cluster_postprocess"]["log"],
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_post_process.R"),
        outdir=cluster.path_tpl["cluster_postprocess"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_low"),
    shell:
        """
        Rscript "{params.script}" \
            --clusters="{input}" \
            --mincells=10 \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule compare_clusters:
    input:
        anndata=cluster.anndata,
        cids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
    output:
        cluster.path_tpl["compare_clusters"]["output"],
    log:
        cluster.path_tpl["compare_clusters"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_compare.py"),
        ncomp="{ncomp}",
        outdir=cluster.path_tpl["compare_clusters"]["dir"],
        reductiontype=cluster.rdim_name,
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input.anndata}" \
            --clusterids="{input.cids}" \
            --ncomp="{params.ncomp}" \
            --outdir="{params.outdir}" \
            --reduced_dims_name="{params.reductiontype}" \
            &> "{log}"
        """


rule clustree:
    input:
        expand(
            cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
            ncomp=["{ncomp}"],
            resolu=cluster.clust_r_lst,
        ),
    output:
        cluster.path_tpl["clustree"]["output"],
    log:
        cluster.path_tpl["clustree"]["log"],
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_clustree.R"),
        res_str=cluster.clust_r_str,
        id_files_str=lambda wc, input: ",".join(input),
        outdir=cluster.path_tpl["clustree"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
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
        neighs=cluster.path_tpl["neighbour_graph"]["output"],
        cids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
        ccolours=cluster.path_tpl["cluster_postprocess"]["outputs"]["ccolors"],
    output:
        cluster.path_tpl["paga"]["outputs"],
    log:
        cluster.path_tpl["paga"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_paga.py"),
        outdir=cluster.path_tpl["paga"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_high"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input.neighs}" \
            --outdir="{params.outdir}" \
            --cluster_ids="{input.cids}" \
            --cluster_colors="{input.ccolours}" \
            &> "{log}"
        """


rule umap:
    input:
        cluster.path_tpl["neighbour_graph"]["output"],
    output:
        cluster.path_tpl["umap"]["output"],
    log:
        cluster.path_tpl["umap"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_umap.py"),
        mindist="{mindist}",
        outdir=cluster.path_tpl["umap"]["dir"],
    threads: 2
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input}" \
            --mindist="{params.mindist}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule plot_rdims_factors:
    input:
        rdims_table=cluster.main_umap_tab,
        metadata=cluster.path_tpl["metadata"]["outputs"]["metatab"],
    output:
        list(cluster.path_tpl["rdims_factors"]["outputs"].values()),
    log:
        cluster.path_tpl["rdims_factors"]["log"],
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_rdims_factor.R"),
        colour_factor_arg=cluster.cfact_arg,
        shape_factor_arg=cluster.sfact_arg,
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=cluster.path_tpl["rdims_factors"]["dir"],
        summary_tex=cluster.path_tpl["rdims_factors"]["outputs"]["summary_tex"],
    resources:
        mem_mb=cluster.get_mem("memory_low"),
    shell:
        """
        Rscript "{params.script}" \
            --table="{input.rdims_table}" \
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
        echo "\\input{{{params.outdir}/UMAP}}" > "{params.summary_tex}"
        """


rule plot_rdims_clusters:
    input:
        rdims_table=cluster.path_tpl["umap"]["output"],
        cluster_ids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
    output:
        png=cluster.path_tpl["rdims_clusters"]["outputs"]["png"],
        tex=cluster.path_tpl["rdims_clusters"]["outputs"]["tex"],
    log:
        cluster.path_tpl["rdims_clusters"]["log"],
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_rdims_factor.R"),
        umap_spec=cluster.umap_spec,
        shape_factor_arg=cluster.sfact_arg,
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=cluster.path_tpl["rdims_clusters"]["dir"],
    resources:
        mem_mb=cluster.get_mem("memory_low"),
    shell:
        """
        Rscript "{params.script}" \
            --method="{params.umap_spec}" \
            --table="{input.rdims_table}" \
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


rule summarise_rdims_clusters:
    input:
        expand(
            rules.plot_rdims_clusters.output.tex,
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            mindist=cluster.mindists_lst,
        ),
    output:
        cluster.path_tpl["rdims_clusters"]["summary_output"],
    shell:
        """
        inputlist=( "{input}" )
        for input in ${{inputlist[@]}}; do
            echo "\\input{{$input}}" >> "{output}"
        done
        """


rule plot_rdims_singler:
    input:
        table=cluster.main_umap_tab,
        labels=cluster.singler_labels_tpl,
    output:
        cluster.path_tpl["rdims_singler"]["output"],
    log:
        cluster.path_tpl["rdims_singler"]["log"],
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_rdims_factor.R",
        reference="{ref}",
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=cluster.path_tpl["rdims_singler"]["dir"],
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


rule plot_singler:
    input:
        metadata=cluster.path_tpl["metadata"]["outputs"]["metatab"],
        scores=cluster.singler_scores_tpl,
        labels=cluster.singler_labels_tpl,
    output:
        cluster.path_tpl["rdims_singler_plot"]["output"],
    log:
        cluster.path_tpl["rdims_singler_plot"]["log"],
    params:
        script=f"{RSCRIPT_DIR}/cluster_singleR_plots.R",
        reference="{ref}",
        outdir=cluster.path_tpl["rdims_singler_plot"]["dir"],
        pdf=cluster.pdf,
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


rule summarise_singler:
    input:
        expand(
            cluster.path_tpl["rdims_singler_plot"]["output"],
            ncomp=["{ncomp}"],
            ref=cluster.singler_ref_lst,
        ),
    output:
        cluster.path_tpl["summarise_singler"]["output"],
    run:
        from utils.cluster import summariseSingleR

        summariseSingleR(
            cluster.path_tpl["summarise_singler"]["dir"],
            cluster.singler_ref_lst,
            output[0],
        )


rule plot_group_numbers:
    input:
        metadata=cluster.path_tpl["metadata"]["outputs"]["metatab"],
        cluster_ids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
    output:
        cluster.path_tpl["group_numbers"]["output"],
    log:
        cluster.path_tpl["group_numbers"]["log"],
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_group_numbers.R"),
        key="{key}",
        options=lambda wc: cluster.populate_options(wc.key),
        outdir=cluster.path_tpl["group_numbers"]["dir"],
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


rule summarise_group_numbers:
    input:
        expand(
            rules.plot_group_numbers.output[0],
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            key=cluster.summary_dict.keys(),
        ),
    output:
        cluster.path_tpl["summarise_group_numbers"]["output"],
    log:
        cluster.path_tpl["summarise_group_numbers"]["log"],
    params:
        outdir=cluster.path_tpl["summarise_group_numbers"]["dir"],
    run:
        from utils.cluster import summariseGroupNumbers

        summariseGroupNumbers(cluster.summary_dict, params.outdir)


# def get_subset_levels(wc):
#     levels_file = cluster.path_tpl["metadata"]["outputs"].get("levels")
#     ck = checkpoints.metadata.get(**wc)

#     if levels_file not in ck.output:
#         return ["all"]
#     else:
#         return pd.read_csv(levels_file, header=None)[0].tolist()


rule cluster_stats:
    input:
        anndata=cluster.anndata,
        cids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
    output:
        stats=cluster.path_tpl["cluster_stats"]["outputs"]["stats"],
        sizes=cluster.path_tpl["cluster_stats"]["outputs"]["sizes"],
    log:
        cluster.path_tpl["cluster_stats"]["log"],
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_stats.py"),
        subset_stat=lambda wc: cluster.subset_stat,
        subset_level="{level}",
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input.anndata}" \
            {params.subset_stat} \
            --subset_level="{params.subset_level}" \
            --clusterids="{input.cids}" \
            --outfile="{output.stats}" \
            &> "{log}"
        """


rule find_markers:
    input:
        anndata=cluster.anndata,
        cids=cluster.path_tpl["cluster_postprocess"]["outputs"]["cids_full"],
        stats=cluster.path_tpl["cluster_stats"]["outputs"]["stats"],
        sizes=cluster.path_tpl["cluster_stats"]["outputs"]["sizes"],
    output:
        cluster.path_tpl["find_markers"]["output"],
    log:
        cluster.path_tpl["find_markers"]["log"],
    params:
        script=f"{PYSCRIPT_DIR}/cluster_markers.py",
        subset_stat=cluster.subset_stat,
        level="{level}",
        cluster="{cluster}",
        markers_test=cluster.test_method,
        markers_pseudocount=cluster.pseudocount,
    resources:
        mem_mb=cluster.get_mem("memory_standard"),
    shell:
        """
        python "{params.script}" \
            --anndata="{input.anndata}" \
            {params.subset_stat} \
            --subset_level="{params.level}" \
            --clusterids="{input.cids}" \
            --cluster="{params.cluster}" \
            --group_means="{input.stats}" \
            --group_sizes="{input.sizes}" \
            --method="{params.markers_test}" \
            --pseudocount="{params.markers_pseudocount}" \
            --outfile="{output}" \
            &> {log}
        """


rule summarise_markers:
    input:
        metadata=cluster.path_tpl["metadata"]["outputs"]["metatab"],
        cids=cluster.path_tpl["scanpy_cluster"]["outputs"]["cids_full"],
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
    resources:
        mem_mb=config["resources"]["mem_mb"],
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


# rule topMarkerHeatmap:
#     input:
#         loom_file=os.path.join(
#             OUT_DIR, "loom.dir", f"{config['plot']['heatmap_matrix']}.loom"
#         ),
#         cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
#         metadata=os.path.join(OUT_DIR, "metadata.dir", "metadata.tsv.gz"),
#         marker_table=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#     output:
#         os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.heatmap.png"),
#     log:
#         os.path.join(MARKERS_DIR("{ncomp}", "{resolu}"), "topMarkerHeatmap.log"),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_top_marker_heatmap.R"),
#         matrix_loc="matrix",
#         scale=(config["plot"]["heatmap_matrix"] == "log1p"),
#         pdf=config["plot"]["pdf"],
#         subgroup=(
#             f"--subgroup={config['plot']['subgroup']}"
#             if config["plot"]["subgroup"]
#             else ""
#         ),
#         outdir=MARKERS_DIR("{ncomp}", "{resolu}"),
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --loom="{input.loom_file}" \
#             --clusterids="{input.cids}" \
#             --metadata="{input.metadata}" \
#             --matrix_loc="{params.matrix_loc}" \
#             --scale="{params.scale}" \
#             --markers="{input.marker_table}" \
#             --pdf="{params.pdf}" \
#             {params.subgroup} \
#             --outdir="{params.outdir}" \
#             &> "{log}"
#         """
# # NOTE: parallelised
# rule dePlots:
#     input:
#         marker_table=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#     output:
#         os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.png"),
#         os.path.join(
#             DE_PLOTS_DIR("{ncomp}", "{resolu}"), "characterise.degenes.{cluster}.tex"
#         ),
#     log:
#         os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.log"),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_de_plots.R"),
#         cluster="{cluster}",
#         outdir=DE_PLOTS_DIR("{ncomp}", "{resolu}"),
#         pdf=config["plot"]["pdf"],
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --degenes="{input.marker_table}" \
#             --cluster="{params.cluster}" \
#             --outdir="{params.outdir}" \
#             --pdf="{params.pdf}" \
#             --plotdirvar=clusterMarkerDEPlotsDir \
#             &> "{log}"
#         """
# def get_clusters_with_marker(ncomp, resolu):
#     markerfile = checkpoints.summariseMarkers.get(ncomp=ncomp, resolu=resolu).output[0]
#     markers = pd.read_csv(markerfile, sep="\t")
#     markers = markers.loc[(markers["p.adj"] < 0.1) & (markers["p.adj"].notna()), :]
#     clusters_with_markers = list(markers["cluster"].unique())
#     assert (
#         "911" not in clusters_with_markers
#     ), "cluster 911 should not used for markerPlots."
#     return clusters_with_markers
# rule summariseDEPlots:
#     input:
#         expand(
#             os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "dePlots.{cluster}.png"),
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
#         ),
#         expand(
#             os.path.join(
#                 DE_PLOTS_DIR("{ncomp}", "{resolu}"),
#                 "characterise.degenes.{cluster}.tex",
#             ),
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
#         ),
#     output:
#         os.path.join(DE_PLOTS_DIR("{ncomp}", "{resolu}"), "summarised_dePlots.sentinel"),
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         touch "{output}"
#         """
# # NOTE: parallelised
# # NOTE: markers in CGAT version line 1244?
# # NOTE: group_opt not implemented
# # NOTE: violinplot not visible when few clusters => height = max(min(nclusters/20 * 5, 10), 3)
# rule markerPlots:
#     input:
#         marker_table=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#         loom_file=os.path.join(LOOM_DIR(), "log1p.loom"),
#         metadata_file=os.path.join(METADATA_DIR(), "metadata.tsv.gz"),
#         cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
#         rdims_table=os.path.join(
#             UMAP_DIR("{ncomp}"), f"umap.{config['plot']['umap_mindist']}.tsv.gz"
#         ),
#     output:
#         rdims=os.path.join(
#             MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.rdims.png"
#         ),
#         violins=os.path.join(
#             MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.violins.png"
#         ),
#         heatmap=os.path.join(
#             MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.heatmap.png"
#         ),
#     log:
#         os.path.join(
#             MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "marker.plots.{cluster}.log"
#         ),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_marker_plots.R"),
#         scaled_loom=(
#             f"--scaled_loom={os.path.join(LOOM_DIR(), 'X.loom')}"
#             if config["plot"]["heatmap_matrix"] == "X"
#             else ""
#         ),
#         rdims_table=os.path.join(
#             UMAP_DIR("{ncomp}"), f"umap.{config['plot']['umap_mindist']}.tsv.gz"
#         ),
#         cluster="{cluster}",
#         outdir=MARKER_PLOTS_DIR("{ncomp}", "{resolu}"),
#         group_opt="",
#         pdf=config["plot"]["pdf"],
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --markers="{input.marker_table}" \
#             --loom="{input.loom_file}" \
#             {params.scaled_loom} \
#             --metadata="{input.metadata_file}" \
#             --clusterids="{input.cids}" \
#             --rdimstable="{input.rdims_table}" \
#             --cluster="{params.cluster}" \
#             --outdir="{params.outdir}" \
#             {params.group_opt} \
#             --pdf="{params.pdf}" \
#             &> "{log}"
#         """
# rule summariseMarkerPlots:
#     input:
#         expand(
#             os.path.join(
#                 MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "cluster.{cluster}.{plot}.png"
#             ),
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
#             plot=["rdims", "violins", "heatmap"],
#         ),
#     output:
#         os.path.join(
#             MARKER_PLOTS_DIR("{ncomp}", "{resolu}"), "summarised_markerPlots.sentinel"
#         ),
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         touch "{output}"
#         """
# # NOTE: cluster_ids not used.
# # NOTE: minfc and minpadj hard coded? => from yaml
# rule plotMarkerNumbers:
#     input:
#         marker_table=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#         cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv.gz"),
#     output:
#         os.path.join(MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"), "deNumbers.png"),
#     log:
#         os.path.join(
#             MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"), "plotMarkerNumbers.log"
#         ),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_plot_marker_numbers.R"),
#         outdir=MARKER_DE_PLOTS_DIR("{ncomp}", "{resolu}"),
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --degenes="{input.marker_table}" \
#             --clusterids="{input.cids}" \
#             --outdir="{params.outdir}" \
#             --minfc=2 \
#             --minpadj=0.05 \
#             --plotdirvar=clusterMarkerDEPlotsDir \
#             &> "{log}"
#         """
# checkpoint plots:
#     input:
#         rules.compareClusters.output,
#         rules.clustTree.output,
#         rules.paga.output,
#         rules.plotRdimsFactors.output,
#         expand(
#             rules.plot_rdims_clusters.output,
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             mindist=str2list(config["plot"]["umap_mindists"]),
#         ),
#         expand(
#             rules.plotGroupNumbers.output,
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             key=config["summaries"].keys(),
#         ),
#         rules.summariseDEPlots.output,
#         rules.summariseMarkerPlots.output,
#     output:
#         os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "plots.sentinel"),
#     shell:
#         """
#         touch "{output}"
#         """
# # NOTE: replaced the dependence of CellHub API by actual annotation files
# rule genesetAnalysis:
#     input:
#         marker_table=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#     output:
#         expand(
#             os.path.join(
#                 GENESETS_DIR("{ncomp}", "{resolu}"),
#                 "genesets.{cluster}.{files}.tsv.gz",
#             ),
#             ncomp=["{ncomp}"],
#             resolu=["{resolu}"],
#             cluster=["{cluster}"],
#             files=[
#                 "GO.BP",
#                 "GO.CC",
#                 "GO.MF",
#                 "KEGG",
#                 "msigdb_biocarta",
#                 "msigdb_reactome",
#             ],
#         ),
#     log:
#         os.path.join(
#             GENESETS_DIR("{ncomp}", "{resolu}"), "geneset.analysis.{cluster}.log"
#         ),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_geneset_analysis.R"),
#         cluster="{cluster}",
#         universe=os.path.join(
#             MARKERS_DIR("{ncomp}", "{resolu}"), "{cluster}.universe.tsv.gz"
#         ),
#         species=config["geneset"]["species"],
#         ensembl=CELLHUB_ANNOT_ENSEMBL,
#         kegg=CELLHUB_ANNOT_KEGG,
#         gmt_names=",".join(config.get("gmt_files", {}).keys()) or "none",
#         gmt_files=",".join(config.get("gmt_files", {}).values()) or "none",
#         adjpthreshold=config["geneset"]["marker_adjpthreshold"],
#         outdir=GENESETS_DIR("{ncomp}", "{resolu}"),
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --markers="{input.marker_table}" \
#             --universe="{params.universe}" \
#             --species="{params.species}" \
#             --annotation="{params.ensembl}" \
#             --kegg_pathways="{params.kegg}" \
#             --gmt_names="{params.gmt_names}" \
#             --gmt_files="{params.gmt_files}" \
#             --cluster="{params.cluster}" \
#             --adjpthreshold="{params.adjpthreshold}" \
#             --direction=positive \
#             --outdir="{params.outdir}" \
#             &> "{log}"
#         for output in "{output}"; do
#             touch $output
#         done
#         """
# # NOTE: possible to have clusters with no gene sets?
# rule summariseGenesetAnalysis:
#     input:
#         lambda wc: [
#             os.path.join(
#                 GENESETS_DIR(wc.ncomp, wc.resolu),
#                 f"genesets.{cluster}.{filetype}.tsv.gz",
#             )
#             for filetype in [
#                 "GO.BP",
#                 "GO.CC",
#                 "GO.MF",
#                 "KEGG",
#                 "msigdb_biocarta",
#                 "msigdb_reactome",
#             ]
#             for cluster in get_clusters_with_marker(wc.ncomp, wc.resolu)
#         ],
#         cids=os.path.join(CLUSTER_DIR("{ncomp}", "{resolu}"), "cluster_ids.tsv"),
#     output:
#         os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.xlsx"),
#         os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.table.tex"),
#         os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets.figure.tex"),
#     log:
#         os.path.join(
#             GENESETS_DIR("{ncomp}", "{resolu}"), "summarise.geneset.analysis.log"
#         ),
#     params:
#         script=os.path.join(RSCRIPT_DIR, "cluster_geneset_summary.R"),
#         genesets_dir=GENESETS_DIR("{ncomp}", "{resolu}"),
#         gmt_names=",".join(config.get("gmt_files", {}).keys()) or "none",
#         show_detailed=config["geneset"]["show_detailed"],
#         min_genes=config["geneset"]["min_fg_genes"],
#         pvalue_threshold=config["geneset"]["pvalue_threshold"],
#         padjust_method=config["geneset"]["padjust_method"],
#         use_adjusted=config["geneset"]["use_adjusted_pvalues"],
#         min_odds_ratio=config["geneset"]["min_odds_ratio"],
#         show_common=config["geneset"]["show_common"],
#         out_prefix=os.path.join(GENESETS_DIR("{ncomp}", "{resolu}"), "cluster.genesets"),
#         pdf=config["plot"]["pdf"],
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         Rscript "{params.script}" \
#             --genesetdir="{params.genesets_dir}" \
#             --gmt_names="{params.gmt_names}" \
#             --show_detailed="{params.show_detailed}" \
#             --clusters="{input.cids}" \
#             --mingenes="{params.min_genes}" \
#             --pvaluethreshold="{params.pvalue_threshold}" \
#             --padjustmethod="{params.padjust_method}" \
#             --useadjusted="{params.use_adjusted}" \
#             --minoddsratio="{params.min_odds_ratio}" \
#             --showcommon="{params.show_common}" \
#             --outprefix="{params.out_prefix}" \
#             --prefix="genesets" \
#             --plotdirvar="clusterGenesetsDir" \
#             --pdf="{params.pdf}" \
#             &> "{log}"
#         """
# rule latexVars:
#     input:
#         rules.taskSummary.output,
#         rules.plots.output,
#         rules.summariseRdimsClusters.output,
#     output:
#         os.path.join(cluster.latex_dir("{ncomp}", "{resolu}"), "report.vars.sty"),
#     params:
#         ncomp="{ncomp}",
#         resolu="{resolu}",
#     run:
#         from utils.cluster import generate_report_vars
#         generate_report_vars(
#             output[0],
#             cluster,
#             CELLHUB_CODE_DIR,
#             params.ncomp,
#             params.resolu,
#         )
# rule summaryReportSource:
#     input:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "report.vars.sty"),
#     output:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "summary.report.tex"),
#     params:
#         ncomp="{ncomp}",
#         resolu="{resolu}",
#     run:
#         from utils.cluster import generate_summary_report
#         generate_summary_report(
#             output[0],
#             cluster,
#             CELLHUB_CODE_DIR,
#             params.ncomp,
#             params.resolu,
#         )
# rule SummaryReport:
#     input:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "summary.report.tex"),
#     output:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "summaryReport.pdf"),
#     log:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "summaryReport.log"),
#     params:
#         run_dir=cluster.out_dir,
#         compilation_dir=os.path.join(
#             LATEX_DIR("{ncomp}", "{resolu}"), "summary.report.dir"
#         ),
#     resources:
#         mem_mb=config["resources"]["memory_standard"],
#     shell:
#         """
#         rm -rf "{params.compilation_dir}"
#         mkdir "{params.compilation_dir}"
#         cd "{params.run_dir}"
#         pdflatex -output-directory="{params.compilation_dir}" \
#             -draftmode \
#             "{input}" \
#             > "{log}"
#         pdflatex -output-directory="{params.compilation_dir}" \
#             "{input}" \
#             > "{log}"
#         mv "{params.compilation_dir}"/summary.report.pdf "{output}"
#         """
# rule markerReportSource:
#     input:
#         marker_table=os.path.join(
#             cluster.markers_dir("{ncomp}", "{resolu}"), "markers.summary.table.tsv.gz"
#         ),
#         latexvars=os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "report.vars.sty"
#         ),
#     output:
#         os.path.join(LATEX_DIR("{ncomp}", "{resolu}"), "marker.report.tex"),
#     params:
#         ncomp="{ncomp}",
#         resolu="{resolu}",
#     run:
#         from utils.cluster import generate_marker_report
#         generate_marker_report(
#             output[0],
#             input.marker_table,
#             input.latexvars,
#             cluster,
#             CELLHUB_CODE_DIR,
#             params.ncomp,
#             params.resolu,
#         )
# rule markerReport:
#     input:
#         os.path.join(cluster.latex_dir("{ncomp}", "{resolu}"), "marker.report.tex"),
#     output:
#         os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "clusterMarkerReport.pdf"
#         ),
#     log:
#         os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "clusterMarkerReport.log"
#         ),
#     params:
#         compilation_dir=os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "marker.report.dir"
#         ),
#     resources:
#         mem_mb=config["resources"]["memory_standard"],
#     shell:
#         """
#         rm -rf "{params.compilation_dir}"
#         mkdir "{params.compilation_dir}"
#         pdflatex -output-directory="{params.compilation_dir}" \
#             -draftmode \
#             "{input}" \
#             > "{log}"
#         pdflatex -output-directory="{params.compilation_dir}" \
#             "{input}" \
#             > "{log}"
#         mv "{params.compilation_dir}"/marker.report.pdf "{output}"
#         """
# rule export:
#     input:
#         expand(
#             os.path.join(cluster.latex_dir("{ncomp}", "{resolu}"), "{rep}"),
#             ncomp="{ncomp}",
#             resolu="{resolu}",
#             rep=["summaryReport.pdf", "clusterMarkerReport.pdf"],
#         ),
#     output:
#         os.path.join(
#             cluster.reports_dir(), "{ncomp}.comps.{resolu}.res", "export.sentinel"
#         ),
#     params:
#         outdir=os.path.join(cluster.reports_dir(), "{ncomp}.comps.{resolu}.res"),
#         # NOTE: between_testfactor no longer exists?
#         # between_xlsx=f"markers.between.{cluster.get_param('markers')}"
#         summary_report=os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "summaryReport.pdf"
#         ),
#         marker_report=os.path.join(
#             cluster.latex_dir("{ncomp}", "{resolu}"), "clusterMarkerReport.pdf"
#         ),
#         markers=os.path.join(
#             cluster.markers_dir("{ncomp}", "{resolu}"), "markers.summary.table.xlsx"
#         ),
#         genesets=os.path.join(
#             cluster.genesets_dir("{ncomp}", "{resolu}"), "cluster.genesets.xlsx"
#         ),
#         # conditions_marker = ...
#         # conditions_genesets = ...
#     shell:
#         """
#         rm -rf "{params.outdir}"
#         mkdir "{params.outdir}"
#         targets=( "{params.summary_report}" "{params.marker_report}" "{params.markers}" "{params.genesets}" )
#         for target_file in ${{targets[@]}}; do
#             if [ -f $target_file ]; then
#                 bname=$(basename $target_file)
#                 ln -s $target_file {params.outdir}/$bname
#             fi
#         done
#         touch {output}
#         """
# def cellxgene_resolutions():
#     if config["cellxgene"]["resolution"] == "all":
#         return RESOLUTION_LST
#     else:
#         return list(config["cellxgene"]["resolution"])
# def cellxgene_resolution_files(resolu_list, ncomp):
#     return [
#         os.path.join(CLUSTER_DIR(ncomp, resolu), "cluster_ids.tsv.gz")
#         for resolu in resolu_list
#     ]
# rule cellxgene:
#     input:
#         cellxgene_resolution_files(cellxgene_resolutions(), "{ncomp}"),
#         anndata=cluster.get_param("anndata"),
#         umap_path=os.path.join(
#             UMAP_DIR("{ncomp}"), f"umap.{config['plot']['umap_mindist']}.tsv.gz"
#         ),
#     output:
#         os.path.join(RDIM_DIR("{ncomp}"), "cellxgene.h5ad"),
#     log:
#         os.path.join(RDIM_DIR("{ncomp}"), "cellxgene.log"),
#     params:
#         script=os.path.join(PYSCRIPT_DIR, "cluster_cellxgene.py"),
#         obs=config["cellxgene"]["obs"],
#         umap_facet_x=config["cellxgene"]["umap_facet_x"],
#         umap_facet_y=config["cellxgene"]["umap_facet_y"],
#         cluster_names=",".join([f"leiden_r{x}" for x in cellxgene_resolutions()]),
#         cluster_paths=",".join(
#             cellxgene_resolution_files(cellxgene_resolutions(), "{ncomp}")
#         ),
#         cluster_split=config["cellxgene"]["cluster_split"],
#     resources:
#         mem_mb=config["resources"]["mem_mb"],
#     shell:
#         """
#         python "{params.script}" \
#             --source_anndata="{input.anndata}" \
#             --obs="{params.obs}" \
#             --umap="{input.umap_path}" \
#             --umap_facet_x="{params.umap_facet_x}" \
#             --umap_facet_y="{params.umap_facet_y}" \
#             --cluster_paths="{params.cluster_paths}" \
#             --cluster_names="{params.cluster_names}" \
#             --cluster_split="{params.cluster_split}" \
#             --adt=None \
#             --outfile="{output}" \
#             &> "{log}"
#         """
