#' ---
#' title: "Plot alluvial diagrams"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/summarise_entropy.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "summarise_entropy.log"
#' ---
#' ---
#' Plot alluvial diagrams to compare two integration tools.
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(Matrix))
stopifnot(require(dplyr))
stopifnot(require(reshape2))
stopifnot(require(ggplot2))
stopifnot(require(tenxutils))
stopifnot(require(ggalluvial))
stopifnot(require(knitr))
stopifnot(require(gridExtra))
stopifnot(require(RColorBrewer))
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
  # Path to cluster ids post-integration for selected tools.
  "samples_plot" = "",
  
  # Name of folders to output rds, pdf and png files.
  "outdir" = ""
)

options(theme_set(theme_bw()))

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
## ############### (i) Read in data  ####################################### ##
## ######################################################################### ##

#' ## Alluvial plots to show relationship between different integration methods 
#' 
#' If only two methods are included, the alluvial and flow plot will look the 
#' same. \
#' The path to the samples compared here are `r paste(opt$samples_plot, collapse=",")`.
#' 

flog.info("Load cluster assignments ...\n")

cluster_assignments = data.frame()
for (s in strsplit(opt$samples_plot, ",")[[1]]){
  clusters <- readRDS(s)
  folder <- strsplit(s, "/")[[1]][2]
  sample_name <- strsplit(folder, "\\.")[[1]][1]
  clusters <- as.data.frame(clusters)
  colnames(clusters) <- sample_name
  clusters$barcode <- rownames(clusters)
  if (nrow(cluster_assignments) == 0){
    cluster_assignments = clusters
  } else {
    cluster_assignments = dplyr::left_join(cluster_assignments, clusters, by="barcode")
  }
}

melted <- melt(cluster_assignments, id.vars = "barcode")

# order clusters
melted$value <- factor(melted$value, levels = c(min(as.numeric(melted$value)): 
                                                max(as.numeric(melted$value))))

## ######################################################################### ##
## ############### (ii) make alluvial plots  ############################### ##
## ######################################################################### ##

#' ### Alluvial flow diagram for selected tools (post-integration)
#' 
#' Here the color is adjusted in each strata (= changes throughout the plot).


flog.info("Plot alluvial diagrams for selected tools (post-integration) ...")

flog.info("Plot as flow diagram (color adjusted in each strata) ...")

gp <- ggplot(melted, aes(x = variable, stratum = value, alluvium = barcode,
                         fill = value), alpha=0.3) 
gp <- gp + geom_flow(stat = "alluvium", lode.guidance = "frontback") 
gp <- gp + geom_stratum(width = 1/10)  + xlab("Integration tool") + ylab("Number of cells")

save_ggplots(file.path(opt$outdir, "Post_integration_alluvial_selected_tools_geom_flow"), 
             gp)

#+ AlluvialFlow_tools, include=TRUE, fig.height=5, fig.cap="", fig.align="center"
gp
#+ include=FALSE
#'


flog.info("Plot as alluvial diagram (same color throughout) ...")

#' ### Alluvial diagram for selected tools (post-integration)
#' 
#' Here the colors stays the same throughout the plot.


gp <- ggplot(melted, aes(x = variable, stratum = value, alluvium = barcode,
                         fill = value), alpha=0.3) 
gp <- gp + geom_alluvium()
gp <- gp + geom_stratum(width = 1/10)  + xlab("Integration tool") + ylab("Number of cells")

save_ggplots(file.path(opt$outdir, "Post_integration_alluvial_selected_tools_geom_alluvial"), 
             gp)

#+ AlluvialStrata_tools, include=TRUE, fig.height=5, fig.cap="", fig.align="center"
gp
#+ include=FALSE
#'

#sessionInfo()

flog.info("Completed.")
