#' ---
#' title: "Make UMAP plots"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/plot_umap.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "plot_umap.log"
#' ---
#' ---
#' Plot UMAP colored for different variables
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
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Path to the post integration UMAP files.
  "coord_file" = "",
  
  # Cluster assignments for the given resolution.
  "clusterids" = NULL,
  
  # Name of the integration tool used.
  "integration_tool" = "harmony",
  
  # Variables to plot on UMAP (from metadata). Grouping variable is automatically added.
  "plot_vars" = "",
  
  # Name of folders to output rds, pdf and png files.
  "outdir" = ""
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

options(theme_set(theme_bw()))

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

input_coord <- read.table(file = gzfile(opt$coord_file), header=TRUE,
                          sep="\t")
tool = opt$integration_tool

## add cluster assignments if required
if (!is.null(opt$clusterids)){
  clusters = as.data.frame(readRDS(opt$clusterids))
  colnames(clusters) = "cluster_postIntegration"
  clusters$barcode = rownames(clusters)
  
  plot_coord = dplyr::left_join(input_coord, clusters, by="barcode")
  
  variables_plot = c(unlist(strsplit(opt$plot_vars,",")), "cluster_postIntegration")
  print(variables_plot)
} else {
  variables_plot = unlist(strsplit(opt$plot_vars,","))
  plot_coord = input_coord
}



## ######################################################################### ##
## ###################### (ii) Make plots ################################## ##
## ######################################################################### ##

#' ## UMAP plots with different variables for integration tool `r opt$integration_tool`
#'    

flog.info("Start plotting ...")

gglist <- list()

for (v in variables_plot) {
  flog.info(paste0("Making plots for variable: ", v))
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + geom_point(alpha=0.5)
    ggsave(file.path(opt$outdir, paste0("umap.", v, ".pdf")), 
           plot = gp, device = cairo_pdf)
    
    gglist[[v]] <- gp
}


#+ umap_variable, include=TRUE, fig.width=6, fig.height=4, fig.cap="", fig.align="center", results="asis"
for (v in variables_plot) {
  cat("\n### UMAP colored by ", v, "\n ")
  print(gglist[[v]])
  cat("\n")
}
#'



flog.info("Completed")