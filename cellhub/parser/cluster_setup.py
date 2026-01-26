import os
import pandas as pd

from parser.parse_args import parse2int, parse_mem, str2list


class ClusterSetup:
    def __init__(self, config):
        # output directory
        self.outdir = config.get("out_dir", "./")
        os.makedirs(self.outdir, exist_ok=True)

        # General information and resource allocations
        self._set_project_info(config)
        self._set_resources(config)
        self._set_tasks(config)

        # Inputs and shared parameters
        self._set_inputs(config)
        self._set_shared_params(config)
        self._set_rdims_params(config["dimension_reduction"])
        self._set_clust_params(config["clustering"])

        # set task-specific parameters
        self._set_preflight()
        self._set_loom(config["plot"])
        self._set_neighbour_graph(config["neighbour_graph"])
        self._set_visualisation_params(config["plot"])
        self._set_umap(config["plot"])
        self._set_group_numbers(config["summaries"])
        self._set_find_markers(config["markers"])
        self._set_geneset_analysis(config["geneset"])
        self._set_cellxgene(config["cellxgene"])

    # <-------------------------- set parameters --------------------------> #
    def _set_project_info(self, config):
        self.projectname = config["projectname"]
        self.author = config["author"]

    def _set_resources(self, config):
        default_resources = {
            "threads": 1,
            "threads_hnsw": 1,
            "mem_low": 4000,
            "mem_std": 8000,
            "mem_high": 16000,
            "time_short": "00:10:00",
            "time_std": "01:00:00",
            "time_long": "12:00:00",
        }
        resource_dict = config.get("resources", {})
        unknown = set(resource_dict) - set(default_resources)
        if unknown:
            raise KeyError(
                f"Unknown resource categories: {sorted(unknown)}. "
                f"Allowed: {sorted(default_resources)}."
            )
        self.resources = {}
        for key, default in default_resources.items():
            value = resource_dict.get(key, default)
            if key.startswith("mem_"):
                value = parse_mem(value)
            self.resources[key] = value

    def _set_tasks(self, config):
        self.task_dict = config["run"]

    def _set_inputs(self, config):
        # input anndata
        self.anndata = config["anndata"]
        # ensembl and kegg annotations and singleR paths
        cellhub_dict = config["cellhub"]
        cellhub_dir = cellhub_dict.get("cellhub_dir")
        if cellhub_dir is not None:
            self.ensembl = os.path.join(
                cellhub_dir, "annotation.dir", "ensembl.to.entrez.tsv.gz"
            )
            self.kegg = os.path.join(cellhub_dir, "annotation.dir", "kegg.pathways.rds")
            if self.task_dict.get("singleR", False):
                self.singler_dir = os.path.join(cellhub_dir, "api", "singleR")
        else:
            self.ensembl = cellhub_dict["ensembl_annotations"]
            self.kegg = cellhub_dict["kegg_pathways"]
            if self.task_dict.get("singleR", False):
                self.singler_dir = cellhub_dict["singleR_dir"]
        # singleR references
        if self.task_dict.get("singleR", False):
            self.singler_labels_tpl = os.path.join(
                self.singler_dir, r"{ref}", "labels.tsv.gz"
            )
            self.singler_scores_tpl = os.path.join(
                self.singler_dir, r"{ref}", "scores.tsv.gz"
            )
            self.singler_ref_lst = os.listdir(self.singler_dir)
        # genesets GMT file
        self.gmt_dict = config.get("gmt_files", None)
        self.gmtname_lst = list(self.gmt_dict.keys()) if self.gmt_dict else []
        self.gmtname_str = ",".join(self.gmtname_lst) if self.gmtname_lst else "none"
        self.gmtfile_lst = list(self.gmt_dict.values()) if self.gmt_dict else []
        self.gmtfile_str = ",".join(self.gmtfile_lst) if self.gmtfile_lst else "none"

    def _set_shared_params(self, config):
        # ----- markers' conserved factor and clustering subset factor
        self.conserved = config["markers"].get("conserved", False)
        if self.conserved:
            self.conserved_arg = "--conserved"
            self.conserved_fact = config["markers"]["conserved_factor"]
            self.subset_stat = "--subset_factor=" + self.conserved_fact
        else:
            self.conserved_arg = ""
            self.conserved_fact = None
            self.subset_stat = ""
        # ----- Data points
        self.pt_shape = config["plot"].get("shape", None)
        self.sfact_arg = (
            "" if self.pt_shape is None else "--shapefactor=" + self.pt_shape
        )
        self.pt_alpha = config["plot"].get("pointalpha", None)
        self.pt_size = config["plot"].get("pointsize", None)
        self.pt_pch = config["plot"].get("pointpch", None)
        # ----- pdf
        self.pdf = config["plot"]["pdf"]

    def _set_rdims_params(self, rdim_dict):
        self.rdim_name = rdim_dict["rdim_name"]
        self.ncomp_str = rdim_dict["n_components"]
        self.ncomp_lst = str2list(self.ncomp_str)
        self.rdims_dir_tpl = os.path.join(self.outdir, r"out.{ncomp}.comp.dir")
        self.max_rdims = max([parse2int(ncomp) for ncomp in self.ncomp_lst])

    def _set_clust_params(self, clust_dict):
        self.clust_algo = clust_dict["algorithm"]
        self.clust_r_str = clust_dict["resolutions"]
        self.clust_r_lst = str2list(self.clust_r_str)
        self.predef_clust_col = clust_dict.get("clust_col", None)
        self.clust_dir_tpl = os.path.join(self.rdims_dir_tpl, r"cluster.{resolu}.dir")

    # <-------------------------- set module output paths --------------------------> #
    def _set_preflight(self):
        self.preflight_geneids = (
            "--gene_ids" if self.task_dict.get("genesets", False) else ""
        )

    def _set_loom(self, plot_dict):
        self.hm_layer = plot_dict["heatmap_matrix"]
        if self.hm_layer == "log1p":
            self.scale_matrix = True
        else:
            self.scale_matrix = False
            self.layers = [self.hm_layer, "log1p"]

    def _set_neighbour_graph(self, neigh_dict):
        self.neigh_method = neigh_dict["method"]
        self.neigh_metric = neigh_dict["metric"]
        self.n_neigh = neigh_dict["n_neighbors"]
        self.hnsw_threads = neigh_dict.get("threads", 1)
        self.hnsw_fullspeed = neigh_dict.get("full_speed", False)

    def _set_visualisation_params(self, plot_dict):
        # ----- Visualisation groups and subgroup
        self.vis_grps_str = plot_dict.get("groups", "cluster")
        self.vis_grps_lst = str2list(self.vis_grps_str)
        if "cluster" not in self.vis_grps_lst:
            self.vis_grps_lst += ["cluster"]
        self.vis_subgrp = plot_dict.get("subgroup", None)
        if self.vis_subgrp is None:
            self.vis_subgrp_arg = ""
        else:
            self.vis_subgrp_arg = f"--subgroup={self.vis_subgrp}"
        # ----- QC vars to visualise
        self.qcvars_str = plot_dict["qcvars"]
        self.qcvars_lst = str2list(self.qcvars_str)
        # ----- Color factors
        self.cfact_lst = self.vis_grps_lst + self.qcvars_lst
        if self.vis_subgrp is not None:
            self.cfact_lst += list(self.vis_subgrp)
        self.cfact_lst = [x for x in set(self.cfact_lst) if x != "cluster"]
        self.cfact_arg = "--colorfactors=" + ",".join(self.cfact_lst)

    def _set_umap(self, plot_dict):
        self.main_mindist = plot_dict["umap_mindist"]
        self.mindists_str = plot_dict["umap_mindists"]
        self.mindists_lst = str2list(self.mindists_str)
        if str(self.main_mindist) not in str2list(self.mindists_lst):
            raise ValueError(
                f"`umap_mindist` should be a value within `umap_mindists` list: {self.mindist} vs {self.mindists_str}."
            )

    def _set_group_numbers(self, summary_dict):
        self.summary_dict = summary_dict

    def _set_find_markers(self, markers_config):
        self.test_method = markers_config["test"]
        self.pseudocount = markers_config["pseudocount"]
        self.min_pct = markers_config["min_pct"]
        self.min_fc = markers_config["min_fc"]
        self.min_padj = markers_config["min_padj"]

    def _set_geneset_analysis(self, geneset_dict):
        self.geneset_names = ["GO.BP", "GO.CC", "GO.MF", "KEGG"] + self.gmtname_lst
        self.species = geneset_dict["species"]
        self.marker_padjthres = geneset_dict["marker_adjpthreshold"]
        self.padj_method = geneset_dict["padjust_method"]
        self.use_padj = geneset_dict["use_adjusted_pvalues"]
        self.pvalthres = geneset_dict["pvalue_threshold"]
        self.min_oddsr = geneset_dict["min_odds_ratio"]
        self.min_fg = geneset_dict["min_fg_genes"]
        self.show_common = geneset_dict["show_common"]
        self.show_detailed = geneset_dict["show_detailed"]

    def _set_cellxgene(self, cellxgene_dict):
        self.cxg_obs = cellxgene_dict["obs"]
        cxg_r = cellxgene_dict["resolution"]
        self.cxg_r = self.clust_r_lst if cxg_r == "all" else str2list(cxg_r)
        self.cxg_facetx = cellxgene_dict["umap_facet_x"]
        self.cxg_facety = cellxgene_dict["umap_facet_y"]
        self.clustsplit = cellxgene_dict["cluster_split"]

    # <-------------------------- utility functions --------------------------> #
    def populate_options(self, summary_key):
        options = []
        for k, v in self.summary_dict[summary_key].items():
            if v == "None" or v == None or v == False or k == "title":
                pass
            elif v == True:
                options.append("--" + k)
            elif k in ["xlab", "ylab"]:
                options.append("--" + k + '="' + str(v) + '"')
            else:
                options.append("--" + k + '="' + str(v) + '"')
        return "\t".join(options)
