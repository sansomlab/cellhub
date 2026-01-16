import os

from parser.parse_args import parse2int, parse_mem, str2list


class ClusterSetup:

    # <-------------------------- load parameters --------------------------> #
    def _load_shared_params(self, config):
        # shared parameters
        self.task_dict = config["run"]
        self.mem_settings = config["resources"]
        if config["markers"].get("conserved", False):
            self.conserved = "--conserved"
            self.conserved_fact = config["markers"]["conserved_factor"]
        else:
            self.conserved = ""
            self.conserved_fact = None
        self.layers = set(str2list(config["plot"]["heatmap_matrix"])) | {"log1p"}

    def _load_cellhub_inputs(self, cellhub_dict):
        cellhub_dir = cellhub_dict.get("cellhub_dir")
        if cellhub_dir is not None:
            self.ensembl = os.path.join(
                cellhub_dir, "annotation.dir", "ensembl.to.entrez.tsv.gz"
            )
            self.kegg = os.path.join(cellhub_dir, "annotation.dir", "kegg.pathways.rds")
            if self.task_dict.get("singleR", False):
                self.singleR_dir = os.path.join(cellhub_dir, "singleR.dir")
        else:
            self.ensembl = cellhub_dict["ensembl_annotations"]
            self.kegg = cellhub_dict["kegg_pathways"]
            if self.task_dict.get("singleR", False):
                self.singleR_dir = cellhub_dict["singleR_dir"]

    def _load_rdims_params(self, rdim_dict):
        self.rdim_name = rdim_dict["rdim_name"]
        self.ncomp_str = rdim_dict["n_components"]
        self.ncomp_lst = str2list(self.ncomp_str)
        self.rdims_dir_tpl = os.path.join(self.outdir, r"out.{ncomp}.comp.dir")
        self.max_rdims = max([parse2int(ncomp) for ncomp in self.ncomp_lst])

    def _load_clust_params(self, clust_dict):
        self.clust_algo = clust_dict["algorithm"]
        self.clust_r_str = clust_dict["resolutions"]
        self.clust_r_lst = str2list(self.clust_r_str)
        self.predef_clust_col = clust_dict.get("clust_col", None)
        self.clust_dir_tpl = os.path.join(self.rdims_dir_tpl, r"cluster.{resolu}.dir")

    # <-------------------------- set module output paths --------------------------> #
    def _set_task_summary(self):
        path_tpl = {
            "log": os.path.join(self.outdir, "task.summary.log"),
            "output": os.path.join(self.outdir, "task.summary.table.tex"),
        }
        self.path_tpl["task_summary"] = path_tpl
        self.target_lst = [path_tpl["output"]]

    def _set_preflight(self):
        self.preflight_geneids = (
            "--gene_ids" if self.task_dict.get("genesets", False) else ""
        )
        path_tpl = {
            "log": os.path.join(self.outdir, "preflight.log"),
        }
        self.path_tpl["preflight"] = path_tpl
        self.target_lst.append(path_tpl["log"])

    def _set_metadata(self):
        metadir = os.path.join(self.outdir, "metadata.dir")
        path_tpl = {
            "dir": metadir,
            "outputs": [os.path.join(metadir, "metadata.tsv.gz")],
            "log": os.path.join(metadir, "metadata.log"),
        }
        if self.conserved:
            path_tpl["outputs"].append(
                os.path.join(metadir, f"{self.conserved_fact}.levels")
            )
        self.path_tpl["metadata"] = path_tpl
        self.target_lst.extend(path_tpl["outputs"])

    def _set_loom(self):
        loomdir = os.path.join(self.outdir, "loom.dir")
        path_tpl = {
            "dir": loomdir,
            "output": os.path.join(loomdir, r"{layer}.loom"),
            "log": os.path.join(loomdir, r"{layer}.log"),
        }
        self.path_tpl["loom"] = path_tpl
        self.target_lst.extend(
            [path_tpl["output"].format(layer=layer) for layer in self.layers]
        )

    def _set_neighbour_graph(self, neigh_dict):
        self.neigh_method = neigh_dict["method"]
        self.neigh_metric = neigh_dict["metric"]
        self.n_neigh = neigh_dict["n_neighbors"]
        self.hnsw_threads = neigh_dict.get("threads", 1)
        self.hnsw_fullspeed = neigh_dict.get("full_speed", False)

        neighdir = os.path.join(self.outdir, "neighbour_graph.dir")
        path_tpl = {
            "dir": neighdir,
            "output": os.path.join(self.rdims_dir_tpl, "neighbour.graph.h5ad"),
            "log": os.path.join(self.rdims_dir_tpl, "neighbour.graph.log"),
        }
        self.path_tpl["neighbour_graph"] = path_tpl
        self.target_lst.extend(
            [path_tpl["output"].format(ncomp=ncomp) for ncomp in self.ncomp_lst]
        )

    def _set_scanpy_cluster(self):
        path_tpl = {
            "dir": self.clust_dir_tpl,
            "output": os.path.join(self.clust_dir_tpl, "scanpy.clusters.tsv.gz"),
            "log": os.path.join(self.clust_dir_tpl, "scanpy.clusters.log"),
        }
        self.path_tpl["scanpy_cluster"] = path_tpl
        self.target_lst.extend(
            [
                path_tpl["output"].format(ncomp=ncomp, resolu=resolu)
                for ncomp in self.ncomp_lst
                for resolu in self.clust_r_lst
            ]
        )

    def _set_cluster_postprocess(self):
        path_tpl = {
            "dir": self.clust_dir_tpl,
            "outputs": {
                "cids_uq": os.path.join(self.clust_dir_tpl, "cluster_ids.tsv"),
                "cids_full": os.path.join(self.clust_dir_tpl, "cluster_ids.tsv.gz"),
                "ccolors": os.path.join(self.clust_dir_tpl, "cluster_colors.tsv"),
                "cccounts": os.path.join(self.clust_dir_tpl, "cluster_cell_counts.tsv"),
            },
            "log": os.path.join(self.clust_dir_tpl, "cluster_postprocess.log"),
        }
        self.path_tpl["cluster_postprocess"] = path_tpl
        self.target_lst.extend(
            [
                outfile.format(ncomp=ncomp, resolu=resolu)
                for ncomp in self.ncomp_lst
                for resolu in self.clust_r_lst
                for outfile in path_tpl["outputs"].values()
            ]
        )

    def _set_compare_clusters(self):
        path_tpl = {
            "dir": self.clust_dir_tpl,
            "output": os.path.join(self.clust_dir_tpl, "cluster.dendrogram.png"),
            "log": os.path.join(self.clust_dir_tpl, "compare_clusters.log"),
        }
        self.path_tpl["compare_clusters"] = path_tpl
        if self.task_dict.get("compare_clusters", False):
            self.target_lst.extend(
                [
                    path_tpl["output"].format(ncomp=ncomp, resolu=resolu)
                    for ncomp in self.ncomp_lst
                    for resolu in self.clust_r_lst
                ]
            )

    def _set_clustree(self):
        path_tpl = {
            "dir": self.rdims_dir_tpl,
            "output": os.path.join(self.rdims_dir_tpl, "clustree.png"),
            "log": os.path.join(self.rdims_dir_tpl, "clustree.log"),
        }
        self.path_tpl["clustree"] = path_tpl
        if self.clust_algo != "predefined":
            self.target_lst.extend(
                [path_tpl["output"].format(ncomp=ncomp) for ncomp in self.ncomp_lst]
            )

    def _set_paga(self):
        path_tpl = {
            "dir": os.path.join(self.clust_dir_tpl, "paga.dir"),
            "outputs": [
                os.path.join(self.clust_dir_tpl, "paga.dir", "draw_graph_fa.png"),
                os.path.join(self.clust_dir_tpl, "paga.dir", "paga.png"),
            ],
            "log": os.path.join(self.clust_dir_tpl, "paga.dir", "paga.log"),
        }
        self.path_tpl["paga"] = path_tpl
        if self.task_dict.get("paga", False):
            self.target_lst.extend(
                [
                    outfile.format(ncomp=ncomp, resolu=resolu)
                    for ncomp in self.ncomp_lst
                    for resolu in self.clust_r_lst
                    for outfile in path_tpl["outputs"]
                ]
            )

    def _set_umap(self, umap_dict):
        self.mindist = umap_dict["umap_mindist"]
        self.mindists_str = umap_dict["umap_mindists"]
        self.mindists_lst = str2list(self.mindists_str)
        if str(self.mindist) not in str2list(self.mindists_lst):
            raise ValueError(
                f"`umap_mindist` should be a value within `umap_mindists` list: {self.mindist} vs {self.mindists_str}."
            )
        umap_dir = os.path.join(self.rdims_dir_tpl, "umap.dir")
        path_tpl = {
            "dir": umap_dir,
            "output": os.path.join(umap_dir, r"umap.{mindist}.tsv.gz"),
            "log": os.path.join(umap_dir, r"umap.{mindist}.log"),
        }
        self.path_tpl["umap"] = path_tpl
        self.target_lst.extend(
            [
                path_tpl["output"].format(ncomp=ncomp, mindist=mdist)
                for ncomp in self.ncomp_lst
                for mdist in self.mindists_lst
            ]
        )

    # <-------------------------- utility functions --------------------------> #
    def print_targets(self):
        print("Target outputs:")
        for target in self.target_lst:
            print(f" - {target}")

    def __init__(self, config):
        # input anndata
        self.anndata = config["anndata"]

        # output directory
        self.outdir = config["out_dir"]
        os.makedirs(self.outdir, exist_ok=True)

        # load parameters
        self._load_shared_params(config)
        self._load_cellhub_inputs(config["cellhub"])
        self._load_rdims_params(config["dimension_reduction"])
        self._load_clust_params(config["clustering"])
        self.path_tpl = dict()

        # set mandatory tasks
        self._set_task_summary()
        self._set_preflight()
        self._set_metadata()
        self._set_loom()
        self._set_neighbour_graph(config["neighbour_graph"])
        self._set_scanpy_cluster()
        self._set_cluster_postprocess()
        self._set_clustree()
        self._set_umap(config["plot"])

        # set optional tasks
        self._set_compare_clusters()
        self._set_paga()

    #     # Summary plots
    #     self.sum_plots = config["summaries"]

    # # Markers
    # tmp = config["markers"]
    # self.test_method = tmp["test"]
    # self.pseudocount = tmp["pseudocount"]
    # self.min_pct = tmp["min_pct"]
    # self.min_fc = tmp["min_fc"]

    #     # Geneset
    #     self.gmt_file_dict = config["gmt_files"]
    #     tmp = config["geneset"]
    #     self.species = tmp["species"]
    #     self.marker_padjthres = tmp["marker_adjpthreshold"]
    #     self.padj_method = tmp["padjust_method"]
    #     self.use_padj = tmp["use_adjusted_pvalues"]
    #     self.pvalthres = tmp["pvalue_threshold"]
    #     self.min_oddsr = tmp["min_odds_ratio"]
    #     self.min_fg = tmp["min_fg_genes"]
    #     self.show_common = tmp["show_common"]
    #     self.show_detailed = tmp["show_detailed"]

    #     # Plots
    #     tmp = config["plot"]
    #     # ----- UMAP

    #     # ----- Visualisation groups and subgroup
    #     self.vis_grps_str = tmp.get("groups", "cluster")
    #     self.vis_grps_lst = str2list(self.vis_grps_str)
    #     if "cluster" not in self.vis_grps_lst:
    #         self.vis_grps_lst += ["cluster"]
    #     self.vis_subgrp = tmp.get("subgroup", None)
    #     # ----- QC vars to visualise
    #     self.qcvars_str = tmp["qcvars"]
    #     self.qcvars_lst = str2list(self.qcvars_str)
    #     # ----- Color factors
    #     self.cfact_lst = self.vis_grps_lst + self.qcvars_lst
    #     if self.vis_subgrp is not None:
    #         self.cfact_lst += list(self.vis_subgrp)
    #     self.cfact_lst = [x for x in set(self.cfact_lst) if x != "cluster"]
    #     # ----- Data points
    #     self.pt_shape = tmp["shape"]
    #     self.pt_alpha = tmp["pointalpha"]
    #     self.pt_size = tmp["pointsize"]
    #     self.pt_pch = tmp["pointpch"]
    #     # ----- Heatmap
    #     self.hm_mat_lst = str2list(tmp["heatmap_matrix"])
    #     self.layers_lst = list(set(["log1p"] + self.hm_mat_lst))
    #     # ----- pdf
    #     self.pdf = tmp["pdf"]

    #     # CellxGene
    #     tmp = config["cellxgene"]
    #     self.cxg_obs = tmp["obs"]
    #     self.cxg_r = tmp["resolution"]
    #     self.cxg_facetx = tmp["umap_facet_x"]
    #     self.cxg_facety = tmp["umap_facet_y"]
    #     self.clustsplit = tmp["cluster_split"]

    #     # <-------------------------- output structure --------------------------> #
    #     self.preflight_paths = {"log": os.path.join(self.outdir, "preflight.log")}
    #     self.task_summary_paths = {
    #         "tex": os.path.join(self.outdir, "task.summary.table.tex")
    #     }
    #     self.metadata_paths = {}

    # <-------------------------- target outputs --------------------------> #
    # def additional_targets(self):
    #     additional_outputs = []
    #     if self.tasks.get("paga", False):
    #         additional_outputs += [
    #             os.path.join(self.paga_paths(ncomp, resolu)["dir"], fig)
    #             for ncomp in self.ncomp_lst
    #             for resolu in self.clust_r_lst
    #             for fig in ["draw_graph_fa.png", "paga.png"]
    #         ]
    #     return additional_outputs

    # if self.tasks.get("compare_clusters", False):
    #     outputs += [
    #         os.path.join(self.cluster_dir(ncomp, resolu), "cluster.dendrogram.png")
    #         for ncomp in self.ncomp_lst
    #         for resolu in self.clust_r_lst
    #     ]
    # if self.task.get("pata", False):
    #     outputs += [
    #         os.path.join(self.paga_dir(ncomp, resolu), fig)
    #         for ncomp in self.ncomp_lst
    #         for resolu in self.clust_r_lst
    #         for fig in ["draw_graph_fa.png", "paga.png"]
    #     ]
    # if self.tasks.get("top_marker_heatmap", False):
    #     outputs += [
    #         os.path.join(
    #             self.markers_dir(ncomp, resolu), "markers.summary.heatmap.png"
    #         )
    #         for ncomp in self.ncomp_lst
    #         for resolu in self.clust_r_lst
    #     ]
    # if self.tasks.get("genesets", False):
    #     outputs += [
    #         os.path.join(self.genesets_dir(ncomp, resolu), "cluster.genesets.xlsx")
    #         for ncomp in self.ncomp_lst
    #         for resolu in self.clust_r_lst
    #     ]

    def rdims_factors_paths(self, ncomp, abspath=True):
        rdims_factors_dir = os.path.join(
            self.rdims_dir(ncomp, abspath), "rdims.visualisation.dir"
        )
        return {
            "dir": rdims_factors_dir,
            "output": os.path.join(rdims_factors_dir, "rdims_factors.tsv"),
            "log": os.path.join(rdims_factors_dir, "rdims_factors.log"),
        }

    def reports_dir(self, abspath=True):
        if abspath:
            return os.path.join(self.outdir, "reports.dir")
        else:
            return "reports.dir"

    def rdims_vis_cluster_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(
            self.cluster_dir(ncomp, resolu, abspath), "rdims.visualisation.dir"
        )

    def group_numbers_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(
            self.cluster_dir(ncomp, resolu, abspath), "group.numbers.dir"
        )

    def stats_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "stats.dir")

    def markers_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "markers.dir")

    def markers_plots_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(
            self.cluster_dir(ncomp, resolu, abspath), "marker.plots.dir"
        )

    def de_plots_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "de.plots.dir")

    def marker_de_plots_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(
            self.cluster_dir(ncomp, resolu, abspath), "marker.de.plots.dir"
        )

    def genesets_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "genesets.dir")

    def latex_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "latex.dir")

    # <-------------------------- memory allocation --------------------------> #
    def get_mem(self, cat_name):
        if cat_name not in self.mem_settings.keys():
            raise KeyError(
                f"Unknown memory category {cat_name}, must be one of {self.mem_settings.keys()}"
            )
        return parse_mem(self.mem_settings[cat_name])


if __name__ == "__main__":
    import yaml

    with open("../yaml/config_cluster.yml", "r") as f:
        config = yaml.safe_load(f)
    setup = ClusterSetup(config)

    assert (
        setup.latex_dir(15, 0.3, False) == "out.15.comp.dir/cluster.0.3.dir/latex.dir"
    ), f"{setup.latex_dir(15, 0.3, False)} != 'out.15.comp.dir/cluster.0.3.dir/latex.dir'"

    assert (
        setup.get_mem("task_summary") == 8000
    ), f"{setup.get_mem('task_summary')} != 8000"

    assert (
        setup.get_mem("top_marker_heatmap") == 16000
    ), f"{setup.get_mem('top_marker_heatmap')} != 16000"

    assert setup.get_mem("paga") == 64000, f"{setup.get_mem('paga')} != 64000"

    assert (
        setup.get_param("plot", "umap_mindists") == "0,0.5,0.7"
    ), f"{setup.get_param('plot', 'umap_mindists')} != 0,0.5,0.7"

    assert setup.get_param("plot", "umap_mindists", True) == [
        "0",
        "0.5",
        "0.7",
    ], f"{setup.get_param('plot', 'umap_mindists', True)} != [0, 0.5, 0.7]"

    print("config_cluster.py passed the test.")
