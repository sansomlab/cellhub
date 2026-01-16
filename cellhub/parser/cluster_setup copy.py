import os

from parser.parse_args import parse_mem, str2list


class ClusterSetup:
    def __init__(self, config):
        # input anndata
        self.anndata = config["anndata"]

        # output directory
        self.outdir = config["outdir"]
        os.makedirs(self.outdir, exist_ok=True)

        # CellHub annotations
        cellhub_dir = config.get("cellhub_dir", None)
        self.ensembl = (
            config["cellhub_ensembl_annotations"]
            if cellhub_dir is None
            else os.path.join(cellhub_dir, "annotation.dir", "ensembl.to.entrez.tsv.gz")
        )
        self.kegg = (
            config["cellhub_kegg_pathways"]
            if cellhub_dir is None
            else os.path.join(cellhub_dir, "annotation.dir", "kegg.pathways.rds")
        )

        # Provided GMT files
        self.gmt_file_dict = config["gmt_files"]

        # Modules to run
        self.tasks = config["run"]

        # Memory categories
        self.mem_settings = config["resources"]

        # Dimension reduction
        rdim_dict = config["dimension_reduction"]
        self.rdim_name = rdim_dict["rdim_name"]
        self.ncomp_str = rdim_dict["n_components"]
        self.ncomp_lst = str2list(self.ncomp_str)

        # Neighbour graph
        neigh_dict = config["neighbor_graph"]
        self.neigh_method = neigh_dict["method"]
        self.neigh_metric = neigh_dict["metric"]
        self.n_neigh = neigh_dict["n_neighbors"]
        self.hnsw_threads = neigh_dict.get("threads", None)
        self.hnsw_fullspeed = neigh_dict.get("full_speed", False)

        # Clustering
        clust_dict = config["clustering"]
        self.clust_algo = clust_dict["algorithm"]
        self.clust_r_str = clust_dict["resolutions"]
        self.clust_r_lst = str2list(self.clust_r_str)
        self.clust_predefined = clust_dict.get("predefined_clusters", False)

        # Summary plots
        self.sum_plots = config["summaries"]

        # Markers
        marker_dict = config["markers"]
        self.test_method = marker_dict["test"]
        self.pseudocount = marker_dict["pseudocount"]
        self.min_pct = marker_dict["min_pct"]
        self.min_fc = marker_dict["min_fc"]
        if marker_dict["conserved"]:
            self.conserved = {
                "factor": marker_dict["conserved_factor"],
                "padj": marker_dict["conserved_padj"],
                "between": marker_dict["conserved_between"],
                "between_factor": marker_dict["conserved_between_factor"],
                "between_padj": marker_dict["conserved_between_padj"],
            }
        else:
            self.conserved = None

        # Geneset
        gset_dict = config["geneset"]
        self.species = gset_dict["species"]
        self.marker_padjthres = gset_dict["marker_adjpthreshold"]
        self.padj_method = gset_dict["padjust_method"]
        self.use_padj = gset_dict["use_adjusted_pvalues"]
        self.pvalthres = gset_dict["pvalue_threshold"]
        self.min_oddsr = gset_dict["min_odds_ratio"]
        self.min_fg = gset_dict["min_fg_genes"]
        self.show_common = gset_dict["show_common"]
        self.show_detailed = gset_dict["show_detailed"]

        # Plots
        plot_dict = config["plot"]
        # ----- UMAP
        self.mindist = plot_dict["umap_mindist"]
        self.mindists_str = plot_dict["umap_mindists"]
        self.mindists_lst = str2list(self.mindists_str)
        if str(self.mindist) not in str2list(self.mindists_lst):
            raise ValueError(
                f"`umap_mindist` should be a value within `umap_mindists` list: {self.mindist} vs {self.mindists_str}."
            )
        # ----- Visualisation groups and subgroup
        self.vis_grps_str = plot_dict["groups"]
        self.vis_grps_lst = str2list(self.vis_grps_str)
        self.vis_subgrps = plot_dict["subgroup"]
        # ----- QC vars to visualise
        self.qcvars_str = plot_dict["qcvars"]
        self.qcvars_lst = str2list(self.qcvars_str)
        # Data points
        self.pt_shape = plot_dict["shape"]
        self.pt_alpha = plot_dict["pointalpha"]
        self.pt_size = plot_dict["pointsize"]
        self.pt_pch = plot_dict["pointpch"]
        # ----- Heatmap
        self.hm_mat = plot_dict["heatmap_matrix"]
        # ----- pdf
        self.pdf = plot_dict["pdf"]

        # CellxGene
        cxg_dict = config["cellxgene"]
        self.cxg_obs = cxg_dict["obs"]
        self.cxg_r = cxg_dict["resolution"]
        self.cxg_facetx = cxg_dict["umap_facet_x"]
        self.cxg_facety = cxg_dict["umap_facet_y"]
        self.clustsplit = cxg_dict["cluster_split"]

    # <-------------------------- output structure --------------------------> #
    def metadata_dir(self, abspath=True):
        if abspath:
            return os.path.join(self.out_dir, "metadata.dir")
        else:
            return "metadata.dir"

    def loom_dir(self, abspath=True):
        if abspath:
            return os.path.join(self.out_dir, "loom.dir")
        else:
            return "loom.dir"

    def comp_dir(self, ncomp, abspath=True):
        if abspath:
            return os.path.join(self.out_dir, f"out.{ncomp}.comp.dir")
        else:
            return f"out.{ncomp}.comp.dir"

    def reports_dir(self, abspath=True):
        if abspath:
            return os.path.join(self.out_dir, "reports.dir")
        else:
            return "reports.dir"

    def rdims_vis_factor_dir(self, ncomp, abspath=True):
        return os.path.join(self.comp_dir(ncomp, abspath), "rdims.visualisation.dir")

    def umap_dir(self, ncomp, abspath=True):
        return os.path.join(self.comp_dir(ncomp, abspath), "umap.dir")

    def cluster_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.comp_dir(ncomp, abspath), f"cluster.{resolu}.dir")

    def rdims_vis_cluster_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(
            self.cluster_dir(ncomp, resolu, abspath), "rdims.visualisation.dir"
        )

    def paga_dir(self, ncomp, resolu, abspath=True):
        return os.path.join(self.cluster_dir(ncomp, resolu, abspath), "paga.dir")

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
        return parse_mem(self.mem_settings[cat_name])

    # <-------------------------- config parameter --------------------------> #
    def get_param(self, param_name1, param_name2=None, parse2list=False, default=None):
        if param_name1 not in self.params:
            raise KeyError(f"unknown field name: {param_name1}.")

        if param_name2 is None:
            return self.params[param_name1]

        if param_name2 not in self.params[param_name1]:
            raise KeyError(
                f"unknown subfield name under {param_name1} field: {param_name2}."
            )

        if not parse2list:
            return self.params[param_name1][param_name2]
        else:
            return str2list(self.params[param_name1][param_name2])


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
