#' ---
#' title: "Plots for harmony integration"
#' output:
#'  html_document:
#'   self_contained: false
#'   toc: true
#'   toc_float: true
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "integrate_data.log"
#' ---
#' ---
#' Perform integration of a single-cell dataset
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(knitr))
stopifnot(require(futile.logger))

# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path,
	       dev.args = list(png = list(type = "cairo")))

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Name of folders to output rds, pdf and png files
  "outdir" = "",
  
  # Path to the seurat begin.rds object; this needs to be
  # normalized and scaled already
  "seurat_obj" = "",
  
  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",
  
  # Path to the file with highly variable features that should be used for
  # the PCA + harmony integration
  "hv_genes" = NULL,
  
  # Harmony or rawdata
  "tool" = "harmony",
  
  # Variable for harmony, 30 used for rawdata
  "nPCs" = 30,  
  
  # (Harmony) Width of soft kmeans clusters, default = 0.1
  "sigma" = 0.1,
  
  # (Harmony) Diversity clustering penalty parameter, default = 2
  "theta" = 2,
  
  # (Harmony) Ridge regression penalty parameter, default = 1
  "lambda" = 1
)

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

flog.info("Running with options: ", opt, capture = TRUE)
flog.info("\n")

## ######################################################################### ##
## ############# (i) Read in info needed for hvg plots ##################### ##
## ######################################################################### ##

if (endsWith(opt$hv_genes, ".gz")) {
  hvg_info <- read.csv(gzfile(opt$hv_genes), header = TRUE)
} else {
  hvg_info <- read.csv(opt$hv_genes, header = TRUE)
}

features <- hvg_info$gene_name[hvg_info$use_integration == TRUE]

flog.info("Finished reading the input hv genes...")

if (opt$tool == "harmony"){
  txt = paste0(opt$nPCs, " principle components are generated and input into harmony." )
} else {
  txt = "30 principle components are generated and used for umap." 
}

#' ## PCA plots on highly variables genes 
#' The genes provided can be found in `r opt$hv_genes`. \
#' The number of highly variable genes provided is `r length(features)` genes.\
#' `r txt`

## ######################################################################### ##
## ############### (ii) Setion for data integration ######################## ##
## ######################################################################### ##

#' ## Data integration using: `r opt$tool`

# unclear if harmony dims.use argument is working, therefore components defined
# above in RunPCA statement and then full 'pca' slot passed to harmony.
# issue here: https://github.com/immunogenomics/harmony/issues/82

if (opt$tool %in% c("harmony")){
  features_used = paste0("The integration was performed using ", 
                         opt$nPCs," PCs for harmony")
} else if (opt$tool == "rawdata") {
  features_used = "No integration was performed"
}

#' `r features_used`.
#'


## ######################################################################### ##
## #################### (iii) UMAP and elbow plots ######################### ##
## ######################################################################### ##

flog.info("Reading UMAP coordinates ...")
n_comp = 30

plot_data <- read.table(gzfile(file.path(opt$outdir, "umap.tsv.gz")),
                   sep="\t", header = TRUE)

##### UMAP PLOTS
gp <- ggplot(plot_data, aes_string(x="UMAP_1", y="UMAP_2", col=opt$split_var))
gp <- gp +  geom_point(alpha=0.5, size=2) + theme_bw()
gp <- gp + guides(col = guide_legend(override.aes = list(alpha=1)))

ggsave(file.path(opt$outdir, "integrated_umap_nofacet.pdf"), 
       plot = gp, device = cairo_pdf,
       units="in")

#' UMAP (`r n_comp` components) colored by variable that was used for integration.
#+ umap, include=TRUE, fig.height=4, fig.width=6
print(gp)

gp <- gp + facet_wrap(as.formula(paste("~", opt$split_var)), ncol = 3)

w <- min(3, length(unique(plot_data[, opt$split_var])))*4
h <- max(1, length(unique(plot_data[, opt$split_var]))/3)*4

ggsave(file.path(opt$outdir, "integrated_umap_facet_splitvar.pdf"), 
       plot = gp, device = cairo_pdf, width = w, height = h,
       units="in")

opts_chunk$set(fig.width= w, fig.height = h)
#' UMAP (`r n_comp` components) split by the variable that was used for integration.
#+ umap_facet, include=TRUE
print(gp)
#+ include=FALSE

if (opt$tool %in% c("harmony")) {
  if (file.exists(file.path(opt$outdir, paste0(opt$tool,"_stdev.tsv.gz")))) {
    df_stdev = read.table(gzfile(file.path(opt$outdir, paste0(opt$tool,"_stdev.tsv.gz"))),
                          sep="\t", header = TRUE)
    txt <- "Elbow plot for harmony stdev shown below."
    elbow_plot <- ggplot(df_stdev, aes(x=comp, y=stdev)) + geom_point()
    elbow_plot <- elbow_plot + theme_bw() + xlab("Component") + ylab("Stdev")
    h=4
    w=4
  } else {
    flog.info("cannot make elbow plot from harmony dimensions because stdev is absent.")
    h=1
    w=6
    elbow_plot <- ggplot() + theme_void()
    txt <- "No standard deviation was present, therefore no elbow plots."
  }
} else if (opt$tool %in% c("rawdata")) {
  txt <- "No standard deviation was present, thus no elbow plots can be created."
  h=0.1
  w=6
  elbow_plot <- ggplot() + theme_void()
} else {
  flog.info("Incorrect tool given.")
}

# text for elbow plots
#' ## Elbow plots
#'
#' `r txt`
#'
opts_chunk$set(fig.width = w, fig.height = h)
#+ elbowplots, include=TRUE
elbow_plot

flog.info("Completed")



