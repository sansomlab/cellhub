#' ---
#' title: "Run iLISI"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "lisi.log"
#' ---
#' ---
#' Run iLISI post-integration
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(Seurat))
stopifnot(require(Matrix))
stopifnot(require(dplyr))
stopifnot(require(reshape2))
stopifnot(require(tenxutils))
stopifnot(require(knitr))
stopifnot(require(viridis))
stopifnot(require(lisi))
stopifnot(require(gridExtra))
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
  
  # Path to the file with umap coordinates.
  "umap_coord" = "",
  
  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",
  
  # Tool used for integration
  "tool" = "seurat_cca"
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

umap_coord <- read.table(opt$umap_coord, sep="\t", header = TRUE)

n_levels <- length(unique(umap_coord[,c(opt$split_var)]))

vis_cols = c("UMAP_1", "UMAP_2")
vis1 = "UMAP_1"
vis2 = "UMAP_2"
vis = "umap"
outfile_name = "iLisi_output_on_UMAP"


## ######################################################################### ##
## ######################### (ii) Run LISI ################################# ##
## ######################################################################### ##

#' ## Assessment of integration using LISI for batches (iLISI)
#' 
#' The R package lisi is part of the harmony tool. 
#' The paper can be found [here](https://www.nature.com/articles/s41592-019-0619-0).
#' The maximum iLISI is the number of levels the variables used for integration has, 
#' e.g. if there were two batches, the maximum iLISI is 2.     
#' 
#' The tool used here is `r opt$tool`, therefore lisi is run on `r vis` coordinates.    
#' 


# Lisi requires a table with 2 dimensional values per cell plus its metadata
lisi_input <- umap_coord[,c(vis_cols, opt$split_var)]

lisi_input <- lisi_input %>% dplyr::sample_frac(1L, FALSE) %>%
              tidyr::gather(key, val, opt$split_var)

lisi_res <- lisi::compute_lisi(umap_coord[,vis_cols], 
                               umap_coord[,c(opt$split_var), drop = FALSE], opt$split_var)
colnames(lisi_res) <- "lisi_value"

lisi_output <- cbind(umap_coord[,c(vis_cols, "barcode", opt$split_var)], lisi_res) 
lisi_plot <- lisi_output %>% dplyr::sample_frac(1L, FALSE) %>% 
                tidyr::gather(key, lisi_value, lisi_value) 


## ######################################################################### ##
## ######################### (ii) Plot LISI output ######################### ##
## ######################################################################### ##

gp <- ggplot(lisi_plot, aes_string(vis1, vis2, color = "lisi_value")) + geom_point(shape = 21) 
gp <- gp + scale_color_viridis_c() + theme_bw()

save_ggplots(file.path(opt$outdir, outfile_name), 
             gp)

#' Results of iLISI plotted on top of `r vis` coordinates. 
#+ lisi_coordinates, include=TRUE, fig.height=4
print(gp)
#+ include=FALSE

gp <- gp + facet_wrap(as.formula(paste("~", opt$split_var)), ncol=3)

w <- min(3, length(unique(lisi_plot[, opt$split_var])))*4
h <- max(1, length(unique(lisi_plot[, opt$split_var]))/3)*4

opts_chunk$set(fig.width= w, fig.height = h)
#' Results of iLISI plotted on top of `r vis` coordinates, split by variable. 
#+ lisi_coordinates_facet, include=TRUE
print(gp)
#+ include=FALSE

gp <- ggplot(lisi_plot, aes(x=lisi_value)) + geom_histogram()
gp <- gp + coord_cartesian(xlim=c(1, n_levels)) + theme_bw()
gp_histo <- gp

save_ggplots(file.path(opt$outdir, "iLisi_histogram"), 
             gp)

gp <- ggplot(lisi_plot, aes(x=lisi_value)) + geom_density()
gp <- gp + coord_cartesian(xlim=c(1, n_levels)) + theme_bw()
gp_density <- gp

save_ggplots(file.path(opt$outdir, "iLisi_density"), 
             gp)


#' Distribution of iLISI results. 
#+ lisi_distribution, include=TRUE, fig.height=3
grid.arrange(gp_density, gp_histo, ncol = 2)
#+ include=FALSE
#'

gp <- ggplot(lisi_plot, aes_string(x="lisi_value", fill = opt$split_var)) + geom_histogram(show.legend = FALSE)
gp <- gp + coord_cartesian(xlim=c(1, n_levels)) + theme_bw()
gp <- gp + facet_wrap(as.formula(paste("~", opt$split_var)), ncol = 3)

w <- min(3, length(unique(lisi_plot[, opt$split_var])))*2.5
h <- max(1, length(unique(lisi_plot[, opt$split_var]))/3)*2.5

#' Distribution of iLISI by `r opt$split_var`.
opts_chunk$set(fig.width= w, fig.height = h)
#+ iLisi_facet, include=TRUE
print(gp)
#+ include=FALSE


colnames(lisi_output)[colnames(lisi_output) == opt$split_var] <- "iLisi_split_var"
write.csv(lisi_output, file.path(opt$outdir, "ilisi_values_split_var.csv"), row.names = FALSE, quote = FALSE)

flog.info("Completed")
