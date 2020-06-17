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
stopifnot(require(gridExtra))

# set chunk options
# opts_chunk$set(echo=FALSE,
#                warning=FALSE,
#                message = FALSE,
#                include = FALSE,
#                fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Path to the post integration UMAP files.
  "coord_file" = "",
  
  # Variables to plot on UMAP (from metadata). Grouping variable is automatically added.
  "plot_vars" = "",
  
  # Name of folders to output rds, pdf and png files.
  "outdir" = "",
  
  #  Path to metadata file.
  "metadata" = NULL
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

#https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

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

## read in metadata
if (!is.null(opt$metadata)){
  metadata <- read.table(gzfile(opt$metadata), sep="\t", header=TRUE)
  plot_coord <- merge(plot_coord, metadata, by="barcode")
}

## ######################################################################### ##
## ###################### (ii) Make plots ################################## ##
## ######################################################################### ##

#' ## UMAP plots with different variables for integration tool `r opt$integration_tool`
#'    

flog.info("Start plotting ...")

#gglist <- list()

for (v in variables_plot) {
  flog.info("Making plots for variable: %s", v)
  nlevels = length(unique(plot_coord[,v]))
  
  if (nlevels > 15 & !is.numeric(plot_coord[,v])){
    flog.info("Make separate legend file for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + geom_point(alpha=0.5, shape=20) + theme_bw()
    
    # extract legend
    legend <- g_legend(gp)
    gp <- gp + theme(legend.position="none")
    # save legend and plot separately
    ggsave(file.path(opt$outdir, paste0("umap.", v, ".legend.pdf")), 
           plot = legend, device = cairo_pdf)
    ggsave(file.path(opt$outdir, paste0("umap.", v, ".pdf")), 
           plot = gp, device = cairo_pdf)
    #gglist[[v]] <- gp
  } else {
    flog.info("Make plot including legend for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + geom_point(alpha=0.5, shape=20) + theme_bw()
    ggsave(file.path(opt$outdir, paste0("umap.", v, ".pdf")), 
           plot = gp, device = cairo_pdf)
    #gglist[[v]] <- gp    
  }
}

flog.info("Make faceted plots for all variables")

for (v in variables_plot) {
  flog.info(paste0("Making plots for variable: %s", v))
  nlevels = length(unique(plot_coord[,v]))
  
  if (nlevels > 9 & !is.numeric(plot_coord[,v])){
    ncols = floor(sqrt(nlevels))
    flog.info("Make separate legend file for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + geom_point(alpha=0.5, shape=20) + theme_bw()
    if(is.numeric(plot_coord[,v])) {
      gp <- gp + scale_color_viridis_c()
    }
    
    # extract legend
    legend <- g_legend(gp)
    gp <- gp + theme(legend.position="none")
    gp <- gp + facet_wrap(as.formula(paste0("~", v)), ncol=ncols)

    if (nlevels > 70){
      gp <- gp + theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
    }
    
    # set height and width for plot (2 for each panel + 2 for legend)
    w <- min(3, ncols)*4 
    h <- max(1, nlevels/ncols)*4
    
    # save legend and plot separately
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v, ".legend.pdf")), 
           plot = legend, device = cairo_pdf)
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v, ".pdf")), 
           plot = gp, device = cairo_pdf)
    #gglist[[v]] <- gp
  } else if (is.numeric(plot_coord[,v])) {
    flog.info("No faceting required as numeric value.")
  } else {
    flog.info("Make plot including legend for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + geom_point(alpha=0.5, shape=20) + theme_bw()
    gp <- gp + facet_wrap(as.formula(paste0("~", v)), ncol=3)

    # set height and width (2 for each panel + 2 for legend)
    w <- min(3, nlevels)*4 + 2
    h <- max(1, nlevels/3)*4
    
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v, ".pdf")), 
           plot = gp, device = cairo_pdf, height = h, width = w)
    #gglist[[v]] <- gp    
  }
}


#+ umap_variable, include=TRUE, fig.width=6, fig.height=4, fig.cap="", fig.align="center", results="asis"
# for (v in variables_plot) {
#   cat("\n### UMAP colored by ", v, "\n ")
#   print(gglist[[v]])
#   cat("\n")
# }
#'



flog.info("Completed")