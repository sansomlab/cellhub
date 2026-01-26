import os
import sys
import pandas as pd


from parser.parse_args import parse2int, str2list, parse2float
from parser.cluster_setup import ClusterSetup

CELLHUB_CODE_DIR = os.path.abspath(os.path.join(workflow.basedir, os.pardir))
sys.path.append(CELLHUB_CODE_DIR)

RSCRIPT_DIR = os.path.join(CELLHUB_CODE_DIR, "R", "scripts")
PYSCRIPT_DIR = os.path.join(CELLHUB_CODE_DIR, "python")
LATEX_DIR = os.path.join(CELLHUB_CODE_DIR, "cellhub", "reports")

cluster = ClusterSetup(config)

OUTDIR = cluster.outdir
RDIMS_DIR_TPL = r"out.{ncomp}.comp.dir"
CLUSTER_DIR_TPL = os.path.join(RDIMS_DIR_TPL, r"cluster.{resolu}.dir")
RELAT_DIR_TPL = {
    "cellhub_code_dir": CELLHUB_CODE_DIR,
    "rdims": RDIMS_DIR_TPL,
    "cluster": CLUSTER_DIR_TPL,
    "metadata": "metadata.dir",
    "loom": "loom.dir",
    "hm_singler": "singleR.dir",
    "reports": "reports.dir",
    "neighbour_graph": os.path.join(RDIMS_DIR_TPL, "neighbor_graph.dir"),
    "rdims_factors": os.path.join(RDIMS_DIR_TPL, "rdims.visualisation.dir"),
    "rdims_singler": os.path.join(RDIMS_DIR_TPL, "singleR.dir"),
    "umap": os.path.join(RDIMS_DIR_TPL, "umap.dir"),
    "rdims_clusters": os.path.join(CLUSTER_DIR_TPL, "rdims.visualisation.dir"),
    "group_numbers": os.path.join(CLUSTER_DIR_TPL, "group.numbers.dir"),
    "stats": os.path.join(CLUSTER_DIR_TPL, "stats.dir"),
    "paga": os.path.join(CLUSTER_DIR_TPL, "paga.dir"),
    "markers": os.path.join(CLUSTER_DIR_TPL, "markers.dir"),
    "de_plots": os.path.join(CLUSTER_DIR_TPL, "de.plots.dir"),
    "marker_plots": os.path.join(CLUSTER_DIR_TPL, "marker.plots.dir"),
    "marker_de_plots": os.path.join(CLUSTER_DIR_TPL, "marker.de.plots.dir"),
    "genesets": os.path.join(CLUSTER_DIR_TPL, "genesets.dir"),
    "latex": os.path.join(CLUSTER_DIR_TPL, "latex.dir"),
}
DIR_TPL = {k: os.path.join(OUTDIR, v) for k, v in RELAT_DIR_TPL.items()}


# <------------------------- Exit Rules -------------------------> #
rule core:
    input:
        preflight=os.path.join(OUTDIR, "preflight.log"),
        # reports=expand(
        #     os.path.join(DIR_TPL["reports"], "{ncomp}.comps.{resolu}.res", "{report}"),
        #     ncomp=cluster.ncomp_lst,
        #     resolu=cluster.clust_r_lst,
        #     report=[
        #         "summaryReport.pdf",
        #         "clusterMarkerReport.pdf",
        #         "markers.summary.table.xlsx",
        #         "cluster.genesets.xlsx",
        #     ],
        # ),
        cellxgene=expand(
            os.path.join(RDIMS_DIR_TPL, "cellxgene.h5ad"), ncomp=cluster.ncomp_lst
        ),


rule optional:
    input:
        compare_clusters=(
            expand(
                os.path.join(CLUSTER_DIR_TPL, "cluster.dendrogram.png"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
            )
            if cluster.task_dict.get("compare_clusters", False)
            else []
        ),
        paga=(
            expand(
                os.path.join(DIR_TPL["paga"], "{outfile}"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
                outfile=[
                    "paga.png",
                    "draw_graph_fa.paga.initialised.png",
                    "paga_init_fa2.tsv.gz",
                    "umap.paga.initialised.png",
                    "umap.paga.init.tsv.gz",
                ],
            )
            if cluster.task_dict.get("paga", False)
            else []
        ),
        markers_summary=(
            expand(
                os.path.join(DIR_TPL["markers"], "markers.summary.table.tsv.gz"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
            )
            if cluster.task_dict.get("characterise_markers", False)
            else []
        ),
        singler_summary=(
            os.path.join(DIR_TPL["hm_singler"], "summary.tex")
            if cluster.task_dict.get("singleR", False)
            else []
        ),
        top_marker_heatmap=(
            expand(
                os.path.join(DIR_TPL["markers"], "markers.summary.heatmap.png"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
            )
            if cluster.task_dict.get("top_marker_heatmap", False)
            else []
        ),
        de_plots=(
            expand(
                os.path.join(DIR_TPL["de_plots"], "characteriseClusterMarkers.tex"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
            )
            if cluster.task_dict.get("de_plots", False)
            else []
        ),
        genesets=(
            expand(
                os.path.join(DIR_TPL["genesets"], "cluster.genesets.figure.tex"),
                ncomp=cluster.ncomp_lst,
                resolu=cluster.clust_r_lst,
            )
            if cluster.task_dict.get("genesets", False)
            else []
        ),


rule full:  # rule full is required to be after core and optional
    input:
        rules.core.input,
        rules.optional.input,


# <------------------------- Process Rules -------------------------> #
rule task_summary:
    output:
        os.path.join(OUTDIR, "task.summary.table.tex"),
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    run:
        tasks = cluster.task_dict.keys()
        runs = cluster.task_dict.values()
        tab = pd.DataFrame.from_dict({"task": tasks, "run": runs})
        tab.to_latex(buf=output[0], index=False, escape=True)


rule preflight:
    input:
        cluster.anndata,
    output:
        os.path.join(OUTDIR, "preflight.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_preflight.py"),
        rdim_name=cluster.rdim_name,
        max_rdims=cluster.max_rdims,
        conserved=cluster.conserved_arg,
        geneids=cluster.preflight_geneids,
        conserved_factor=cluster.conserved_fact,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_short"],
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
        metatab=[os.path.join(DIR_TPL["metadata"], "metadata.tsv.gz")],
        levels=(
            [os.path.join(DIR_TPL["metadata"], f"{cluster.conserved_fact}.levels")]
            if cluster.conserved
            else []
        ),
    log:
        os.path.join(DIR_TPL["metadata"], "metadata.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_metadata.py"),
        conserved=cluster.conserved_arg,
        conserved_factor=cluster.conserved_fact,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input}" \
            {params.conserved} \
            --conserved_factor="{params.conserved_factor}" \
            --outfile="{output[0]}" \
            &> "{log}"
        """


def get_conserved_levels(wc):
    ck = checkpoints.metadata.get(**wc)
    if cluster.conserved:
        levels_path = os.path.join(
            DIR_TPL["metadata"], f"{cluster.conserved_fact}.levels"
        )
        return pd.read_csv(levels_path, header=None)[0].tolist()
    else:
        return ["all"]


rule loom:
    input:
        cluster.anndata,
    output:
        os.path.join(DIR_TPL["loom"], "{layer}.loom"),
    log:
        os.path.join(DIR_TPL["loom"], "loom.{layer}.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_loom.py"),
        layers="{layer}",
        outdir=DIR_TPL["loom"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        os.path.join(DIR_TPL["neighbour_graph"], "neighbour_graph.h5ad"),
    log:
        os.path.join(DIR_TPL["neighbour_graph"], "neighbour_graph.log"),
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
        threads=cluster.resources["threads_hnsw"],
        fullspeedmode="--fullspeed" if cluster.hnsw_fullspeed else "",
    threads: cluster.resources["threads_hnsw"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        os.path.join(CLUSTER_DIR_TPL, "scanpy.clusters.tsv.gz"),
    log:
        os.path.join(CLUSTER_DIR_TPL, "scanpy.clusters.log"),
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
        outdir=CLUSTER_DIR_TPL,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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


checkpoint cluster_postprocess:
    input:
        rules.scanpy_cluster.output,
    output:
        cids_uq=os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv"),
        cids_full=os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv.gz"),
        ccolors=os.path.join(CLUSTER_DIR_TPL, "cluster_colors.tsv"),
        cccounts=os.path.join(CLUSTER_DIR_TPL, "cluster_cell_counts.tsv"),
    log:
        os.path.join(CLUSTER_DIR_TPL, "cluster_postprocess.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_post_process.R"),
        outdir=CLUSTER_DIR_TPL,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --clusters="{input}" \
            --mincells=10 \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


def get_valid_clusters(ncomp, resolu):
    ck = checkpoints.cluster_postprocess.get(ncomp=ncomp, resolu=resolu)
    cid_file = ck.output["cids_uq"]
    return [x for x in pd.read_csv(cid_file, header=None)[0] if x != 911]


rule compare_clusters:
    input:
        anndata=cluster.anndata,
        cids=os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv.gz"),
    output:
        os.path.join(CLUSTER_DIR_TPL, "cluster.dendrogram.png"),
    log:
        os.path.join(CLUSTER_DIR_TPL, "compare_clusters.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_compare.py"),
        ncomp="{ncomp}",
        outdir=CLUSTER_DIR_TPL,
        reductiontype=cluster.rdim_name,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
            os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv.gz"),
            ncomp=["{ncomp}"],
            resolu=cluster.clust_r_lst,
        ),
    output:
        os.path.join(RDIMS_DIR_TPL, "clustree.png"),
    log:
        os.path.join(RDIMS_DIR_TPL, "clustree.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_clustree.R"),
        res_str=cluster.clust_r_str,
        id_files_str=lambda wc, input: ",".join(input),
        outdir=RDIMS_DIR_TPL,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        neighs=os.path.join(DIR_TPL["neighbour_graph"], "neighbour_graph.h5ad"),
        cids=os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv.gz"),
        ccolours=os.path.join(CLUSTER_DIR_TPL, "cluster_colors.tsv"),
    output:
        expand(
            os.path.join(DIR_TPL["paga"], "{outfile}"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            outfile=[
                "paga.png",
                "draw_graph_fa.paga.initialised.png",
                "paga_init_fa2.tsv.gz",
                "umap.paga.initialised.png",
                "umap.paga.init.tsv.gz",
            ],
        ),
    log:
        os.path.join(DIR_TPL["paga"], "paga.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_paga.py"),
        outdir=DIR_TPL["paga"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_high"],
        time=cluster.resources["time_long"],
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
        os.path.join(DIR_TPL["neighbour_graph"], "neighbour_graph.h5ad"),
    output:
        os.path.join(DIR_TPL["umap"], r"umap.{mindist}.tsv.gz"),
    log:
        os.path.join(DIR_TPL["umap"], r"umap.{mindist}.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_umap.py"),
        mindist="{mindist}",
        outdir=DIR_TPL["umap"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        rdims_table=os.path.join(DIR_TPL["umap"], f"umap.{cluster.main_mindist}.tsv.gz"),
        metadata=rules.metadata.output.metatab[0],
    output:
        expand(
            os.path.join(DIR_TPL["rdims_factors"], "{outfile}"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            outfile=[f"UMAP.{fact}.png" for fact in cluster.cfact_lst],
        ),
        os.path.join(DIR_TPL["rdims_factors"], "UMAP.tex"),
        summary_tex=os.path.join(DIR_TPL["rdims_factors"], "plot.rdims.factor.tex"),
    log:
        os.path.join(DIR_TPL["rdims_factors"], "plot.rdims.factor.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_rdims_factor.R"),
        colour_factor_arg=cluster.cfact_arg,
        shape_factor_arg=cluster.sfact_arg,
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=DIR_TPL["rdims_factors"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
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
        echo "\\input{{{params.outdir}/UMAP}}" > "{output.summary_tex}"
        """


rule plot_rdims_clusters:
    input:
        rdims_table=rules.umap.output,
        cluster_ids=rules.cluster_postprocess.output.cids_full,
    output:
        png=os.path.join(
            DIR_TPL["rdims_clusters"], "umap.mindist_{mindist}.cluster_id.png"
        ),
        tex=os.path.join(DIR_TPL["rdims_clusters"], "umap.mindist_{mindist}.tex"),
    log:
        os.path.join(DIR_TPL["rdims_clusters"], "plot.rdims.cluster.{mindist}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_rdims_factor.R"),
        umap_spec="umap.mindist_{mindist}",
        shape_factor_arg=cluster.sfact_arg,
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=DIR_TPL["rdims_clusters"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
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
        os.path.join(DIR_TPL["rdims_clusters"], "plot.rdims.factor.tex"),
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    shell:
        """
        inputlist=( "{input}" )
        for input in ${{inputlist[@]}}; do
            echo "\\input{{$input}}" >> "{output}"
        done
        """


rule plot_rdims_singler:
    input:
        table=rules.plot_rdims_factors.input.rdims_table,
        labels=cluster.singler_labels_tpl,
    output:
        os.path.join(DIR_TPL["rdims_singler"], "UMAP.{ref}.pruned.labels.png"),
    log:
        os.path.join(DIR_TPL["rdims_singler"], "rdims.plots.{ref}.log"),
    params:
        script=f"{RSCRIPT_DIR}/cluster_plot_rdims_factor.R",
        reference="{ref}",
        pointsize=cluster.pt_size,
        pointalpha=cluster.pt_alpha,
        pointpch=cluster.pt_pch,
        pdf=cluster.pdf,
        outdir=DIR_TPL["rdims_singler"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
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
        metadata=rules.metadata.output.metatab,
        labels=cluster.singler_labels_tpl,
        scores=cluster.singler_scores_tpl,
    output:
        os.path.join(DIR_TPL["hm_singler"], r"{ref}.heatmap.png"),
    log:
        os.path.join(DIR_TPL["hm_singler"], r"singleR.plots.{ref}.log"),
    params:
        script=f"{RSCRIPT_DIR}/cluster_singleR_plots.R",
        reference="{ref}",
        outdir=DIR_TPL["hm_singler"],
        pdf=cluster.pdf,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        umaps=expand(
            rules.plot_rdims_singler.output,
            ncomp=cluster.ncomp_lst,
            ref=cluster.singler_ref_lst,
        ),
        heatmaps=expand(
            rules.plot_singler.output,
            ref=cluster.singler_ref_lst,
        ),
    output:
        # NOTE: changed location from rdims_singler to hm_singler
        os.path.join(DIR_TPL["hm_singler"], "summary.tex"),
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    run:
        from utils.cluster import summariseSingleR

        summariseSingleR(DIR_TPL["hm_singler"], cluster.singler_ref_lst, output[0])


rule plot_group_numbers:
    input:
        metadata=rules.metadata.output.metatab,
        cluster_ids=rules.cluster_postprocess.output.cids_full,
    output:
        os.path.join(DIR_TPL["group_numbers"], r"{key}.data.tsv.gz"),
    log:
        os.path.join(DIR_TPL["group_numbers"], r"plot.group.numbers.{key}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_group_numbers.R"),
        key="{key}",
        options=lambda wc: cluster.populate_options(wc.key),
        outdir=DIR_TPL["group_numbers"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
            rules.plot_group_numbers.output,
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            key=cluster.summary_dict.keys(),
        ),
    output:
        os.path.join(DIR_TPL["group_numbers"], "number.plots.tex"),
    log:
        os.path.join(DIR_TPL["group_numbers"], "summarise.group.numbers.log"),
    params:
        outdir=DIR_TPL["group_numbers"],
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    run:
        from utils.cluster import summariseGroupNumbers

        summariseGroupNumbers(cluster.summary_dict, params.outdir)


rule cluster_stats:
    input:
        anndata=cluster.anndata,
        cids=rules.cluster_postprocess.output.cids_full,
    output:
        stats=os.path.join(DIR_TPL["stats"], "{level}.stats.tsv.gz"),
        sizes=os.path.join(DIR_TPL["stats"], "{level}.sizes.tsv.gz"),
    log:
        os.path.join(DIR_TPL["stats"], "{level}.stats.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_stats.py"),
        subset_stat=cluster.subset_stat,
        subset_level="{level}",
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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
        cids=rules.cluster_postprocess.output.cids_full,
        stats=rules.cluster_stats.output.stats,
        sizes=rules.cluster_stats.output.sizes,
    output:
        os.path.join(DIR_TPL["markers"], "{cluster}.{level}.markers.tsv.gz"),
    log:
        os.path.join(DIR_TPL["markers"], "{cluster}.{level}.markers.log"),
    params:
        script=f"{PYSCRIPT_DIR}/cluster_markers.py",
        subset_stat=cluster.subset_stat,
        level="{level}",
        cluster="{cluster}",
        markers_test=cluster.test_method,
        markers_pseudocount=cluster.pseudocount,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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


def get_valid_cluster_marker_files(wc):
    return [
        os.path.join(DIR_TPL["markers"], "{cluster}.{level}.markers.tsv.gz").format(
            ncomp=wc.ncomp, resolu=wc.resolu, level=level, cluster=cluster
        )
        for level in get_conserved_levels(wc)
        for cluster in get_valid_clusters(wc.ncomp, wc.resolu)
    ]


checkpoint summarise_markers:
    input:
        metadata=rules.metadata.output.metatab,
        cids=rules.cluster_postprocess.output.cids_full,
        cluster_files=lambda wc: get_valid_cluster_marker_files(wc),
    output:
        markers_tsv=os.path.join(DIR_TPL["markers"], "markers.summary.table.tsv.gz"),
        markers_xlsx=os.path.join(DIR_TPL["markers"], "markers.summary.table.xlsx"),
        stats_tsv=os.path.join(DIR_TPL["markers"], "markers.summary.stats.tsv"),
    log:
        os.path.join(DIR_TPL["markers"], "markers_summary.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_summarise_markers.R"),
        markers_str=lambda wc: ",".join(get_valid_cluster_marker_files(wc)),
        min_pct=cluster.min_pct,
        min_fc=cluster.min_fc,
        outdir=DIR_TPL["markers"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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


rule top_marker_heatmap:
    input:
        loom_file=os.path.join(DIR_TPL["loom"], f"{cluster.hm_layer}.loom"),
        cids=rules.cluster_postprocess.output.cids_full,
        metadata=rules.metadata.output.metatab,
        marker_table=rules.summarise_markers.output.markers_tsv,
    output:
        os.path.join(DIR_TPL["markers"], "markers.summary.heatmap.png"),
    log:
        os.path.join(DIR_TPL["markers"], "topMarkerHeatmap.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_top_marker_heatmap.R"),
        matrix_loc="matrix",
        scale=cluster.scale_matrix,
        pdf=cluster.pdf,
        subgroup=cluster.vis_subgrp_arg,
        outdir=DIR_TPL["markers"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
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


rule de_plots:
    input:
        mkgtab=rules.summarise_markers.output.markers_tsv,
    output:
        png=os.path.join(DIR_TPL["de_plots"], r"dePlots.{cluster}.png"),
        tex=os.path.join(DIR_TPL["de_plots"], r"characterise.degenes.{cluster}.tex"),
    log:
        os.path.join(DIR_TPL["de_plots"], r"dePlots.{cluster}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_de_plots.R"),
        cluster="{cluster}",
        outdir=DIR_TPL["de_plots"],
        pdf=cluster.pdf,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --degenes="{input.mkgtab}" \
            --cluster="{params.cluster}" \
            --outdir="{params.outdir}" \
            --pdf="{params.pdf}" \
            --plotdirvar=clusterMarkerDEPlotsDir \
            &> "{log}"
        """


def get_clusters_with_marker(ncomp, resolu):
    markerfile = checkpoints.summarise_markers.get(ncomp=ncomp, resolu=resolu).output[0]
    markers = pd.read_csv(markerfile, sep="\t")
    clusters_with_markers = (
        markers.loc[(markers["p.adj"] < 0.1) & (markers["p.adj"].notna()), "cluster"]
        .unique()
        .tolist()
    )
    assert (
        "911" not in clusters_with_markers
    ), "cluster 911 should not used for markerPlots."
    return clusters_with_markers


rule summarise_de_plots:
    input:
        expand(
            rules.de_plots.output.tex,
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
        ),
    output:
        os.path.join(DIR_TPL["de_plots"], "characteriseClusterMarkers.tex"),
    shell:
        """
        for texpath in {input}; do
            texfile=$(basename "$texpath")
            echo "\\input{{\\clusterMarkerDEPlotsDir/$texfile}}" >> "{output}"
        done
        """


rule marker_plots:
    input:
        mkgtab=rules.summarise_markers.output.markers_tsv,
        loom_file=os.path.join(DIR_TPL["loom"], "log1p.loom"),
        metatab=rules.metadata.output.metatab,
        cids=rules.cluster_postprocess.output.cids_full,
        rdims_tab=os.path.join(DIR_TPL["umap"], f"umap.{cluster.main_mindist}.tsv.gz"),
    output:
        rdims=os.path.join(DIR_TPL["marker_plots"], r"cluster.{cluster}.rdims.png"),
        violins=os.path.join(DIR_TPL["marker_plots"], r"cluster.{cluster}.violins.png"),
        heatmap=os.path.join(DIR_TPL["marker_plots"], r"cluster.{cluster}.heatmap.png"),
    log:
        os.path.join(DIR_TPL["marker_plots"], r"marker.plots.{cluster}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_marker_plots.R"),
        scaled_loom=(
            f"--scaled_loom={os.path.join(DIR_TPL['loom'], 'X.loom')}"
            if cluster.hm_layer == "X"
            else ""
        ),
        rdims_table=os.path.join(DIR_TPL["umap"], f"umap.{cluster.main_mindist}.tsv.gz"),
        cluster="{cluster}",
        outdir=DIR_TPL["marker_plots"],
        group_opt=(f"--group={cluster.vis_subgrp}" if cluster.vis_subgrp else ""),
        pdf=cluster.pdf,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --markers="{input.mkgtab}" \
            --loom="{input.loom_file}" \
            {params.scaled_loom} \
            --metadata="{input.metatab}" \
            --clusterids="{input.cids}" \
            --rdimstable="{input.rdims_tab}" \
            --cluster="{params.cluster}" \
            --outdir="{params.outdir}" \
            {params.group_opt} \
            --pdf="{params.pdf}" \
            &> "{log}"
        """


rule plot_marker_numbers:
    input:
        mkgtab=rules.summarise_markers.output.markers_tsv,
    output:
        os.path.join(DIR_TPL["marker_de_plots"], "deNumbers.png"),
    log:
        os.path.join(DIR_TPL["marker_de_plots"], "plotMarkerNumbers.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_plot_marker_numbers.R"),
        outdir=DIR_TPL["marker_de_plots"],
        minfc=cluster.min_fc,
        minpadj=cluster.min_padj,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --degenes="{input.mkgtab}" \
            --outdir="{params.outdir}" \
            --minfc={params.minfc} \
            --minpadj={params.minpadj} \
            --plotdirvar=clusterMarkerDEPlotsDir \
            &> "{log}"
        """


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


# NOTE: replaced the dependence of CellHub API by actual annotation files
rule geneset_analysis:
    input:
        mkgtab=rules.summarise_markers.output.markers_tsv,
    output:
        expand(
            os.path.join(DIR_TPL["genesets"], "genesets.{cluster}.{files}.tsv.gz"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            cluster=["{cluster}"],
            files=cluster.geneset_names,
        ),
    log:
        os.path.join(DIR_TPL["genesets"], "geneset.analysis.{cluster}.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_geneset_analysis.R"),
        cluster="{cluster}",
        universe=os.path.join(DIR_TPL["markers"], r"{cluster}.universe.tsv.gz"),
        species=cluster.species,
        ensembl=cluster.ensembl,
        kegg=cluster.kegg,
        gmt_names=cluster.gmtname_str,
        gmt_files=cluster.gmtfile_str,
        adjpthreshold=cluster.marker_padjthres,
        outdir=DIR_TPL["genesets"],
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --markers="{input.mkgtab}" \
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
        """


rule summarise_geneset_analysis:
    input:
        lambda wc: [
            os.path.join(DIR_TPL["genesets"], f"genesets.{cluster}.{filetype}.tsv.gz")
            for filetype in cluster.geneset_names
            for cluster in get_clusters_with_marker(wc.ncomp, wc.resolu)
        ],
        cids=rules.cluster_postprocess.output.cids_uq,
    output:
        os.path.join(DIR_TPL["genesets"], "cluster.genesets.xlsx"),
        os.path.join(DIR_TPL["genesets"], "cluster.genesets.table.tex"),
        os.path.join(DIR_TPL["genesets"], "cluster.genesets.figure.tex"),
    log:
        os.path.join(DIR_TPL["genesets"], "summarise.geneset.analysis.log"),
    params:
        script=os.path.join(RSCRIPT_DIR, "cluster_geneset_summary.R"),
        genesets_dir=DIR_TPL["genesets"],
        gmt_names=cluster.gmtname_str,
        show_detailed=cluster.show_detailed,
        min_genes=cluster.min_fg,
        pvalue_threshold=cluster.pvalthres,
        padjust_method=cluster.padj_method,
        use_adjusted=cluster.use_padj,
        min_odds_ratio=cluster.min_oddsr,
        show_common=cluster.show_common,
        out_prefix=os.path.join(DIR_TPL["genesets"], "cluster.genesets"),
        pdf=cluster.pdf,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_std"],
    shell:
        """
        Rscript "{params.script}" \
            --genesetdir="{params.genesets_dir}" \
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


rule latex_vars:
    input:
        task_summary=rules.task_summary.output,
        clustree=rules.clustree.output,
        rdims_factors=rules.plot_rdims_factors.output,
        rdims_clusters=rules.summarise_rdims_clusters.output,
        group_numbers=rules.summarise_group_numbers.output,
        compare_clusters=(
            rules.compare_clusters.output
            if cluster.task_dict.get("compare_clusters")
            else []
        ),
        singler_summary=(
            os.path.join(DIR_TPL["hm_singler"], "summary.tex")
            if cluster.task_dict.get("singleR", False)
            else []
        ),
        paga=(rules.paga.output if cluster.task_dict.get("paga") else []),
        marker_plots=(
            expand(
                rules.marker_plots.output,
                ncomp=["{ncomp}"],
                resolu=["{resolu}"],
                cluster=lambda wc: get_clusters_with_marker(wc.ncomp, wc.resolu),
            )
        ),
        markers_summary=(
            os.path.join(DIR_TPL["markers"], "markers.summary.table.tsv.gz")
            if cluster.task_dict.get("characterise_markers", False)
            else []
        ),
        de_plots=(
            rules.summarise_de_plots.output
            if cluster.task_dict.get("de_plots")
            else []
        ),
        genesets=(
            os.path.join(DIR_TPL["genesets"], "cluster.genesets.figure.tex")
            if cluster.task_dict.get("genesets", False)
            else []
        ),
    output:
        os.path.join(DIR_TPL["latex"], "report.vars.sty"),
    params:
        ncomp="{ncomp}",
        resolu="{resolu}",
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    run:
        from utils.cluster import generate_report_vars

        generate_report_vars(
            output[0], RELAT_DIR_TPL, cluster, params.ncomp, params.resolu
        )


rule summary_report_source:
    input:
        rules.latex_vars.output,
    output:
        os.path.join(DIR_TPL["latex"], "summary.report.tex"),
    params:
        ncomp="{ncomp}",
        resolu="{resolu}",
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    run:
        from utils.cluster import generate_summary_report

        generate_summary_report(
            output[0], cluster, DIR_TPL, params.ncomp, params.resolu
        )


rule summary_report:
    input:
        rules.summary_report_source.output,
    output:
        os.path.join(DIR_TPL["latex"], "summaryReport.pdf"),
    log:
        os.path.join(DIR_TPL["latex"], "summaryReport.log"),
    params:
        run_dir=cluster.outdir,
        compilation_dir=os.path.join(DIR_TPL["latex"], "summary.report.dir"),
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    shell:
        """
        rm -rf "{params.compilation_dir}"
        mkdir "{params.compilation_dir}"
        cd "{params.run_dir}"
        pdflatex -output-directory="{params.compilation_dir}" \
            -draftmode \
            "{input}" \
            > "{log}"
        pdflatex -output-directory="{params.compilation_dir}" \
            "{input}" \
            > "{log}"
        mv "{params.compilation_dir}"/summary.report.pdf "{output}"
        """


rule marker_report_source:
    input:
        marker_table=rules.summarise_markers.output.markers_tsv,
        latexvars=rules.latex_vars.output,
    output:
        os.path.join(DIR_TPL["latex"], "marker.report.tex"),
    params:
        ncomp="{ncomp}",
        resolu="{resolu}",
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    run:
        from utils.cluster import generate_marker_report

        generate_marker_report(
            output[0],
            input.marker_table,
            input.latexvars,
            DIR_TPL,
            cluster,
            params.ncomp,
            params.resolu,
        )


rule marker_report:
    input:
        rules.marker_report_source.output,
    output:
        os.path.join(DIR_TPL["latex"], "clusterMarkerReport.pdf"),
    log:
        os.path.join(DIR_TPL["latex"], "clusterMarkerReport.log"),
    params:
        compilation_dir=os.path.join(DIR_TPL["latex"], "marker.report.dir"),
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    shell:
        """
        rm -rf "{params.compilation_dir}"
        mkdir "{params.compilation_dir}"
        pdflatex -output-directory="{params.compilation_dir}" \
            -draftmode \
            "{input}" \
            > "{log}"
        pdflatex -output-directory="{params.compilation_dir}" \
            "{input}" \
            > "{log}"
        mv "{params.compilation_dir}"/marker.report.pdf "{output}"
        """


rule export:
    input:
        expand(
            os.path.join(DIR_TPL["latex"], "{rep}"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            rep=["summaryReport.pdf", "clusterMarkerReport.pdf"],
        ),
    output:
        expand(
            os.path.join(DIR_TPL["reports"], "{ncomp}.comps.{resolu}.res", "{report}"),
            ncomp=["{ncomp}"],
            resolu=["{resolu}"],
            report=[
                "summaryReport.pdf",
                "clusterMarkerReport.pdf",
                "markers.summary.table.xlsx",
                "cluster.genesets.xlsx",
            ],
        ),
    params:
        outdir=os.path.join(DIR_TPL["reports"], "{ncomp}.comps.{resolu}.res"),
        summary_report=os.path.join(DIR_TPL["latex"], "summaryReport.pdf"),
        marker_report=os.path.join(DIR_TPL["latex"], "clusterMarkerReport.pdf"),
        markers=os.path.join(DIR_TPL["markers"], "markers.summary.table.xlsx"),
        genesets=os.path.join(DIR_TPL["genesets"], "cluster.genesets.xlsx"),
    threads: 1
    resources:
        mem_mb=cluster.resources["mem_low"],
        time=cluster.resources["time_short"],
    shell:
        """
        rm -rf "{params.outdir}"
        mkdir "{params.outdir}"
        targets=( "{params.summary_report}" "{params.marker_report}" "{params.markers}" "{params.genesets}" )
        for target_file in ${{targets[@]}}; do
            if [ -f $target_file ]; then
                bname=$(basename $target_file)
                ln -s $target_file {params.outdir}/$bname
            fi
        done
        """


def cellxgene_resolution_files(resolu_list, ncomp):
    return [
        os.path.join(CLUSTER_DIR(ncomp, resolu), "cluster_ids.tsv.gz")
        for resolu in resolu_list
    ]


rule cellxgene:
    input:
        resolu_files=expand(
            os.path.join(CLUSTER_DIR_TPL, "cluster_ids.tsv.gz"),
            ncomp=["{ncomp}"],
            resolu=cluster.cxg_r,
        ),
        anndata=cluster.anndata,
        umap_path=os.path.join(DIR_TPL["umap"], f"umap.{cluster.main_mindist}.tsv.gz"),
    output:
        os.path.join(RDIMS_DIR_TPL, "cellxgene.h5ad"),
    log:
        os.path.join(RDIMS_DIR_TPL, "cellxgene.log"),
    params:
        script=os.path.join(PYSCRIPT_DIR, "cluster_cellxgene.py"),
        obs=cluster.cxg_obs,
        umap_facet_x=cluster.cxg_facetx,
        umap_facet_y=cluster.cxg_facety,
        cluster_names=",".join([f"leiden_r{x}" for x in cluster.cxg_r]),
        cluster_paths=lambda wc, input: ",".join(input.resolu_files),
        cluster_split=cluster.clustsplit,
    threads: cluster.resources["threads"]
    resources:
        mem_mb=cluster.resources["mem_std"],
        time=cluster.resources["time_std"],
    shell:
        """
        python "{params.script}" \
            --source_anndata="{input.anndata}" \
            --obs="{params.obs}" \
            --umap="{input.umap_path}" \
            --umap_facet_x="{params.umap_facet_x}" \
            --umap_facet_y="{params.umap_facet_y}" \
            --cluster_paths="{params.cluster_paths}" \
            --cluster_names="{params.cluster_names}" \
            --cluster_split="{params.cluster_split}" \
            --adt=None \
            --outfile="{output}" \
            &> "{log}"
        """
