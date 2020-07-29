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
stopifnot(require(dplyr))
stopifnot(require(ggrepel))

# set chunk options
# opts_chunk$set(echo=FALSE,
#                warning=FALSE,
#                message = FALSE,
#                include = FALSE,
#                fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------

# these parameters are passed from Rscript run
option_list <- list(
  make_option(
    c("--task_yml"),
    dest = "task_yml",
    help="Path to yml file"
  ),
  make_option(
    c("--log_filename"),
    dest = "log_filename",
    help="Path to log file"
  ))
params <- parse_args(OptionParser(option_list=option_list))

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))
flog.info("Running with parameters: ", params, capture = TRUE)


# The script expects the following paramters from yml:
default_options <- list(
  # Path to the post integration UMAP files.
  "coord_file" = "",

  # Variables to plot on UMAP (from metadata). Grouping variable is automatically added.
  "plot_vars" = "",

  # If set in pipeline.yml, then path to cluster assignments table is passed to here.
  # Otherwise stays NULL.
  "plot_clusters" = NULL,

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

flog.info("Running with options:", opt, capture = TRUE)
flog.info("\n")

#https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

gp_nolabels <- theme(axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

gp_choose_pointsize <- function(size){
  if(size > 50000){
    g <- geom_point(alpha=0.5, size=1, shape='.')
  } else if (size > 5000) {
    g <- geom_point(alpha=0.5, size=0.5, shape=20)
  } else {
    g <- geom_point(alpha=0.7, size=1, shape=20)
  }
  return(g)
}

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

input_coord <- read.table(file = gzfile(opt$coord_file), header=TRUE,
                          sep="\t")
tool = opt$integration_tool

## add cluster assignments if required
if (!is.null(opt$plot_clusters)){
  clusters = read.table(gzfile(opt$plot_clusters), header=TRUE, sep="\t")
  colnames(clusters) = c("barcode", "cluster_id")

  plot_coord = merge(input_coord, clusters, by="barcode", all.x = TRUE)
  flog.info("Cluster assignments are available for %s cells.",
            length(which(is.na(plot_coord$cluster_id) == FALSE)))
  flog.info("Cluster assignments not available for %s cells.",
            length(which(is.na(plot_coord$cluster_id) == TRUE)))
  plot_coord$cluster_id = factor(plot_coord$cluster_id,
                                 levels = c(min(plot_coord$cluster_id,
                                                na.rm = TRUE):max(plot_coord$cluster_id,
                                                                  na.rm = TRUE),
                                            NA))
  variables_plot = c(unlist(strsplit(opt$plot_vars,",")), "cluster_id")
} else {
  variables_plot = unlist(strsplit(opt$plot_vars,","))
  plot_coord = input_coord
}

## read in metadata
if (!is.null(opt$metadata)){
  metadata <- read.table(gzfile(opt$metadata), sep="\t", header=TRUE, as.is=TRUE)
  plot_coord <- merge(plot_coord, metadata, by="barcode")
}

flog.info("Number of cells in this sample: %s", nrow(input_coord))

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

  if (nlevels > 15 & !is.numeric(plot_coord[,v]) || nrow(input_coord) > 10000){
    flog.info("Make separate legend file for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + gp_choose_pointsize(nrow(input_coord)) + theme_bw()
    if (!is.numeric(plot_coord[,v])){
      if (nlevels > 9){
        ncols_legend = ceiling(nlevels/9)
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10, shape=16), ncol=ncols_legend))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=10, shape=16)))
      } else {
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10, shape=16)))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=10, shape=16)))
      }
    }
    if(is.numeric(plot_coord[,v])) {
      gp <- gp + scale_color_viridis_c()
    }
    if(v == "cluster_id"){
      # make labels df with mean coordinates per cluster
      labels <- plot_coord %>%
        group_by(cluster_id) %>%
        summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      gp <- gp + geom_text_repel(data = labels, aes(x=UMAP_1,y=UMAP_2, label = cluster_id,
                                                    color=1), color='black', show.legend = FALSE)
    }
    # extract legend
    legend <- g_legend(gp)
    gp <- gp + theme(legend.position="none") + gp_nolabels
    # save legend and plot separately
    v_out = gsub("_", "-", v)
    ggsave(file.path(opt$outdir, paste0("umap.", v_out, ".legend.png")),
           plot = legend, type = "cairo")
    ggsave(file.path(opt$outdir, paste0("umap.", v_out, ".png")),
           plot = gp, type = "cairo")
    #gglist[[v]] <- gp
  } else {
    flog.info("Make plot including legend for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + gp_choose_pointsize(nrow(input_coord))
    gp <- gp + theme_bw() + gp_nolabels
    if (!is.numeric(plot_coord[,v])){
      if (nlevels > 9){
        ncols_legend = ceiling(nlevels/9)
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10), ncol=ncols_legend))
      } else {
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10)))
      }
    }
    if(is.numeric(plot_coord[,v])) {
      gp <- gp + scale_color_viridis_c()
    }
    if(v == "cluster_id"){
      # make labels df with mean coordinates per cluster
      labels <- plot_coord %>%
        group_by(cluster_id) %>%
        summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      gp <- gp + geom_text_repel(data = labels, aes(x=UMAP_1,y=UMAP_2, label = cluster_id,
                                                    color=1), color='black', show.legend = FALSE)
    }
    # extract legend
    legend <- g_legend(gp)
    # save legend and plot separately
    v_out = gsub("_", "-", v)
    ggsave(file.path(opt$outdir, paste0("umap.", v_out, ".legend.png")),
           plot = legend, type = "cairo")
    ggsave(file.path(opt$outdir, paste0("umap.", v_out, ".png")),
           plot = gp, type = "cairo")
    #gglist[[v]] <- gp
  }
}

flog.info("Make faceted plots for all variables")

for (v in variables_plot) {
  flog.info(paste0("Making plots for variable: %s", v))
  nlevels = length(unique(plot_coord[,v]))

  if (nlevels > 9 & !is.numeric(plot_coord[,v]) || nrow(input_coord) > 50000 & !is.numeric(plot_coord[,v])){
    ncols = floor(sqrt(nlevels))
    flog.info("Make separate legend file for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + gp_choose_pointsize(nrow(input_coord)) + theme_bw()
    # split legend into columns
    if (nlevels > 9){
      ncols_legend = ceiling(nlevels/9)
      gp <- gp + guides(colour = guide_legend(override.aes = list(size=10), ncol=ncols_legend))
    } else {
      gp <- gp + guides(colour = guide_legend(override.aes = list(size=10)))
    }
    # extract legend
    legend <- g_legend(gp)
    gp <- gp + theme(legend.position="none") + gp_nolabels
    gp <- gp + facet_wrap(as.formula(paste0("~", v)), ncol=ncols)

    if (nlevels > 50){
      gp <- gp + theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
    }

    # set height and width for plot (2 for each panel + 2 for legend)
    w <- min(3, ncols)*4
    h <- max(1, nlevels/ncols)*4

    # save legend and plot separately
    v_out = gsub("_", "-", v)
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v_out, ".legend.png")),
           plot = legend, type = "cairo")
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v_out, ".png")),
           plot = gp, type = "cairo")
    #gglist[[v]] <- gp
  } else if (is.numeric(plot_coord[,v])) {
    flog.info("No faceting required as numeric value.")
  } else {
    flog.info("Make plot including legend for: %s", v)
    gp <- ggplot(plot_coord, aes_string(x="UMAP_1", y="UMAP_2", color = v))
    gp <- gp + gp_choose_pointsize(nrow(input_coord))
    gp <- gp + theme_bw() + gp_nolabels
    if (!is.numeric(plot_coord[,v])){
      if (nlevels > 9){
        ncols_legend = ceiling(nlevels/9)
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10), ncol=ncols_legend))
      } else {
        gp <- gp + guides(colour = guide_legend(override.aes = list(size=10)))
      }
    }
    gp <- gp + facet_wrap(as.formula(paste0("~", v)), ncol=3)

    # set height and width (2 for each panel + 2 for legend)
    w <- min(3, nlevels)*3 + 2
    h <- max(1, nlevels/3)*3

    v_out = gsub("_", "-", v)
    ggsave(file.path(opt$outdir, paste0("umap.facet.", v_out, ".png")),
           plot = gp, type = "cairo", height = h, width = w)
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