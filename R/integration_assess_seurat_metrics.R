#' ---
#' title: "Assess the integration using metrics from the Seurat package"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "assess_integration_seurat_metrics.log"
#' ---
#' ---
#' Run MixingMetric and LocalStructure from Seurat package to assess the 
#' integration
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(viridis))
stopifnot(require(reshape2))
stopifnot(require(knitr))
stopifnot(require(futile.logger))


# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Name of folders to output rds, pdf and png files
  "outdir" = "",
  
  # Path to the seurat object post-integration
  "seurat_obj" = "",
  
  # Tool used for integration
  "integration_tool" = "seurat_cca",

  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",
  
  # Correct assay to use
  "assay" = "integrated",
  
  # Correct dimension reduction slot to use
  "slot" = "pca",
  
  # Maximum size of local neighborhood to compute in 
  # 'MixingMetric' function of Seurat (in publication script = 300)
  "assessment_max_k" = 300,
  
  # Path to UMAP to visualise MM
  "umap_coord" = "",
  
  # Number of neighbors for LocalStruct function (default = 100)
  "assessment_neighbors" = 100
)

# here the yaml can also be read in from the default location in code 
# directory
options <- read_yaml(params$task_yml)

# Update the default options
if(!is.null(options)) {
 opt <- utils::modifyList(default_options, options)
} else{
 opt <- default_options
}

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))

flog.info("Running with options:", opt, capture = TRUE)
flog.info("\n")

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

s_integrated <- readRDS(file = opt$seurat_obj)

flog.info("Default assay of the object is %s", DefaultAssay(s_integrated))

# write out info about assay and slot that were selected
flog.info("Chosen assay is %s", opt$assay)
flog.info("Chosen slot is %s", opt$slot)

# change the default assay if required as MixingMetric function does not have
# an assay argument
if (DefaultAssay(s_integrated) != opt$assay){
  flog.info("Reset the DefaultAssay because it did not match the 
            required assay.")
  DefaultAssay(s_integrated) <- opt$assay
}

## ######################################################################### ##
## ################ (ii) Run Mixing Metric ################################# ##
## ######################################################################### ##

umap_coord <- read.table(gzfile(opt$umap_coord), sep="\t", header = TRUE)


flog.info("Running Mixing metric from Seurat package...")

#' ## Seurat Metrics: Mixing Metric
#' 
#' Data was integrated using: `r opt$integration_tool`. 
#' The assay used for this assessment is `r opt$assay` and the slot is `r opt$slot`.     
#' The higher this value, the better the batches/variables are mixed.


mm = opt$assessment_max_k - MixingMetric(s_integrated, grouping.var = opt$split_var,
                                dims = 1:2, max.k = opt$assessment_max_k,
                                reduction = opt$slot)

s_integrated@meta.data$barcode = rownames(s_integrated[[]])
metrics = data.frame(MixingMetric = mm, barcode = s_integrated[[]]$barcode)

umap_plot = dplyr::left_join(metrics, umap_coord[,c("barcode", "UMAP_1", "UMAP_2")],
                             by = "barcode")
gp <- ggplot(data=umap_plot, aes(x=UMAP_1, y=UMAP_2, 
                                 color=MixingMetric)) + geom_point(size=0.5) 
gp <- gp + theme_bw() + scale_color_viridis_c() 

ggsave(file.path(opt$outdir, "UMAP_colored_Mixing_metric.pdf"), 
       plot = gp, device = cairo_pdf)

#' UMAP colored by Mixing Metric.
#+ plot_umap_mixingmetric, include=TRUE, fig.height=4
print(gp)
#+ include=FALSE

metrics = dplyr::left_join(metrics, s_integrated[[]], by="barcode")
## plots of mixing metric by grouping variable
metrics %>% group_by_at(opt$split_var) %>% summarise(mean(MixingMetric))

gp <- ggplot(data=metrics, aes_string(x=opt$split_var, y="MixingMetric", 
                                      fill=opt$split_var)) 
gp <- gp + geom_violin(scale = "width", show.legend = FALSE)  + theme_bw()
gp <- gp + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave(file.path(opt$outdir, "Mixing_metric_by_grouping_var.pdf"), 
       plot = gp, device = cairo_pdf)

#' Distribution of Mixing Metric colored by the variable that was used for integration.
#+ plot_mixingmetric, include=TRUE, fig.height=4
print(gp)
#+ include=FALSE


## ######################################################################### ##
## ############## (iii) Run Local Structure ################################ ##
## ######################################################################### ##

#' ## Seurat Metrics: Local Structure
#' 
#' Data was integrated using: `r opt$integration_tool`.  
#' The assay used for this assessment is `r opt$assay` and the slot is `r opt$slot`.     
#' The higher this value, the better the local structure is preserved.  

flog.info("Running Local Structure metric from Seurat package...")

if (opt$integration_tool == "rawdata"){
  ## local structure cannot be run for rawdata
  flog.info("Local Structure cannot be assessed for unintegrated data ... ")
  metrics_ls = data.frame(sample = unique(s_integrated[[]][,opt$split_var]),
                          meanLocalStructure = rep(0, length(unique(s_integrated[[]][,opt$split_var]))))
  colnames(metrics_ls) = c(opt$split_var, "meanLocalStructure")
  # make an empty plot for markdown code
  gp <- ggplot() + theme_void()
  w=0.2
  h=0.2
} else {
  flog.info("Assess Local Structure ...")
  #In Seurat script, they set the default assay to RNA but here leave it
  # slot where the integration is located!
  # two opposing scripts from Satija group (only two occurences of LocalStruct): https://github.com/satijalab/Integration2019/blob/e5821bd242fa0a46eb6fd37764275737512032a4/analysis_code/integration/integration_metrics.R
  # and https://github.com/satijalab/Integration2019/blob/e5821bd242fa0a46eb6fd37764275737512032a4/figure_code/mca_figure.R
  #DefaultAssay(s_integrated) <- "RNA"
  ls = LocalStruct(s_integrated, grouping.var = opt$split_var, neighbors = opt$assessment_neighbors,
                   reduction = opt$slot, reduced.dims = 1:15, orig.dims = 1:15)
  metrics_ls = melt(lapply(ls, mean))
  colnames(metrics_ls) = c("meanLocalStructure", opt$split_var)

  plot_metrics = melt(ls)

  gp <- ggplot(data=plot_metrics, aes(x=value, fill=L1)) 
  gp <- gp + geom_histogram(show.legend = FALSE) + theme_bw()
  gp <- gp + facet_wrap(~L1, ncol = 3)
  ggsave(file.path(opt$outdir, "Localstructure_by_grouping_var.pdf"), 
         plot = gp, device = cairo_pdf, width = 9)
  w <- min(3, length(unique(plot_metrics$L1)))*2.5
  h <- max(1, length(unique(plot_metrics$L1))/3)*2.5
}



#' Distribution of Local Structure metric colored by the variable that 
#' was used for integration. Empty plot indicates that it cannot be assessed (rawdata).
opts_chunk$set(fig.width= w, fig.height = h)
#+ plot_localstructure, include=TRUE
print(gp)
#+ include=FALSE
#' 

flog.info("Finished plotting ...")

# make a summary table for this sample
summary_metrics = metrics %>% group_by_at(opt$split_var) %>% summarise(mean_mixing_metric = mean(MixingMetric))
summary_metrics = dplyr::left_join(summary_metrics, metrics_ls, by=opt$split_var)
summary_metrics$method = opt$integration_tool

# write out the summary table
write.csv(summary_metrics, file = file.path(opt$outdir, paste0("metrics.csv")))

# sink(file=params$log_filename, append = TRUE, type = c("output", "message"))
# sessionInfo()
# sink()
flog.info("Completed")
