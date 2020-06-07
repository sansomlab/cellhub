#' ---
#' title: "Run log-normalization and assessment of highly variable genes"
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
#' Perform normalization of a single-cell dataset
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
stopifnot(require(tidyr))
stopifnot(require(dplyr))
stopifnot(require(Matrix))
stopifnot(require(reshape2))
stopifnot(require(future))
stopifnot(require(gridExtra))
stopifnot(require(cowplot))
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

  # Path to the seurat begin.rds object; this needs to be
  # normalized and scaled already (pre_integrated.rds)
  "seurat_obj" = "",

  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",

  # (Harmony) Whether to perform normalisation per batch or on merged
  # object. By default is set to 'merged' according to current tutorial.
  # merged | perbatch
  "merge_normalisation" = "merged",

  # Variables to regress, example: percent.mito
  "regress_latentvars" = NULL,

  # Whether to regress cell cycle, options: none, difference, all
  "regress_cellcycle" = "none",

  # The genelists to annotate variable genes with
  "vargenes_dir" = NULL,
  
  # Number of variable features to use find (e.g. 3000)
  "ngenes" = 2000
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
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

flog.info("Reading RDS file...  ")
s.full <- readRDS(opt$seurat_obj)
flog.info("Default assay of Seurat object: %s", DefaultAssay(s.full))

# Set variables to regress
if ( ! is.null(opt$regress_latentvars)){
  if(grepl(",", opt$regress_latentvars)){
    vreg <- unlist(strsplit(opt$regress_latentvars, split=","))
  } else {
    vreg <- opt$regress_latentvars
  }
} else {
  vreg = NULL
  flog.info("No variable regression...")
}


## make plot pre normalisation if cell cycle is included.

if ( identical(opt$regress_cellcycle, "none") ) {
  flog.info("No plots produced")
  txt <- "No cell cycle correction performed."
  gp_pre <- ggplot() + theme_void()
} else {
  txt <- paste0("PCA plot before and after cell cycle correction with option: ",
                opt$regress_cellcycle)
  gp_pre <- readRDS(file.path(opt$outdir, "gp_cellcyclePre"))
}

if (opt$merge_normalisation %in% c("merged")){
  txt_merged <- "Normalization was performed on the merged object."
} else {
  txt_merged <- paste0("Normalization was performed per batch (using batch variable: ", opt$split_var, ".")
}


## add cell cycle to regression variables
if ( identical(opt$regress_cellcycle, "none") ) {
  flog.info("Data will be scaled without correcting for cell cycle")
} else {
  if ( identical(opt$regress_cellcycle, "all") ){
    flog.info("Cell cycle correction for S and G2M scores will be applied")
    vreg <- c(vreg, "S.Score", "G2M.Score")
  } else if ( identical(opt$regress_cellcycle, "difference") ) {
    flog.info("Cell cycle correction for the difference between G2M and S phase scores will be applied")
    vreg <- c(vreg, "CC.Difference")
  } else {
    stop("Cell cycle regression type not recognised")
  }
}

## ######################################################################### ##
## ########## (ii) Section for normalisation and hvg ####################### ##
## ######################################################################### ##

#' ## Data normalization
#'
#' The data was normalized using log-normalization.
#' `r txt_merged`
#' The following variables were regressed: `r paste(vreg, collapse=",")` .


## ######################################################################### ##
## ###################### (iii) assess hv genes ############################ ##
## ######################################################################### ##

hvg_info <- read.csv(gzfile(file.path(opt$outdir, 
                                      "hv_genes_info.csv.gz")))
gglist <- readRDS(file.path(opt$outdir, "gglist_hvg.rds"))
flog.info("Finished reading of hv genes table and loading plots ...")

if (opt$merge_normalisation %in% c("merged")) {
  max_iter = 1
  loop = FALSE
} else { 
  max_iter = length(unique(hvg_info$sample))
  loop = TRUE
}
  

if(loop) {
  nr <- ceiling(max_iter/3)
  gp <- plot_grid(plotlist = gglist, 
                  nrow=nr)
} else {
  gp <- gglist
}

w <- min(3, max_iter)*4 + 1
h <- max(1, max_iter/3)*4


ggsave(file.path(opt$outdir, "var_genes_plot.pdf"), 
       plot = gp, device = cairo_pdf,width=w,
       height=h, units="in")

#' ## Assessment of highly-variable genes
#' Number of variable genes selected per sample: `r opt$ngenes` genes.
opts_chunk$set(fig.width= w, fig.height = h)
#+ hvg_genes, include=TRUE
gp

flog.info("Finished plots for hv genes ...")


# add info from genelists that were specified
if (!is.null(opt$vargenes_dir)){
  names_genesets = c()
  for (f in list.files(opt$vargenes_dir, pattern = "*.tsv*")){
    if (endsWith(f, ".gz")){
      list_name = substr(f, 1, nchar(f)-nchar(".tsv.gz"))
    } else {
      list_name = substr(f, 1, nchar(f)-nchar(".tsv"))
    }
  names_genesets = c(names_genesets, list_name)
  }
  flog.info("Make plots of fraction of hvg within each given geneset.")
  df_plot = hvg_info[, c("sample", names_genesets)] %>% pivot_longer(all_of(names_genesets)) %>%
            group_by(sample, name) %>% summarise(n_set = n(),
                                                 fraction_geneset = length(which(value == TRUE))/n_set) 
  nsamples = length(unique(df_plot$sample))
  if(opt$merge_normalisation %in% c("perbatch")){ 
    df_plot$sample = factor(as.character(df_plot$sample), levels = names(gglist))
    }
  gp <- ggplot(df_plot, aes(x=name, y=fraction_geneset)) + geom_bar(stat='identity', fill = "grey")
  gp <- gp + facet_wrap(~sample, ncol = min(nsamples, 3)) + ylab("Fraction of geneset \n within highly variable genes")
  gp <- gp + xlab("Geneset") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                                                  strip.background = element_blank())
  #gp <- gp + coord_cartesian(ylim = c(0,1))
  
  w <- min(3, max_iter)*4
  h <- max(1, max_iter/3)*4 + 2
  
  ggsave(file.path(opt$outdir, "hvg_fraction_genelist.pdf"), 
         plot = gp, device = cairo_pdf, width=w,
         height=h, units="in")
  
  txt_hvg = "Fraction of highly variable genes within different genesets defined in pipeline.yml.\n The x axis labels reflect the genesets provided, while each facets shows a sample/batch."
} else {
  flog.info("No genelists provided")
  gp <- ggplot() + theme_void()
  h = 0.1
  txt_hvg = "No genesets to assess in highly variable genes provided."
}

#' `r txt_hvg`
opts_chunk$set(fig.width= w, fig.height = h)
#+ hvg_genes_lists, include=TRUE
gp



#### check cell cycle regression

if ( identical(opt$regress_cellcycle, "none") ) {
  flog.info("No cell cycle plots produced post correction...")
  gp_post <- ggplot() + theme_void()
  h = 0.1
  w = 4
} else { 
  if (opt$merge_normalisation %in% c("merged")){
    gp_post <- DimPlot(object = s.full, group.by="Phase",
                       reduction = "pcacellcyclePost")
    w=4
  } else {
    gglist <- readRDS(file.path(opt$outdir, "gglist_cellcyclePost.rds"))
    nr <- ceiling(length(gglist)/3)
    gp_post <- plot_grid(plotlist = gglist,  nrow=nr)
    w=12
  }
  h=4
}

# text for cell cycle correct here (to include both plots)
#' ## Cell cycle correction
#'
#' `r txt`
#'
opts_chunk$set(fig.width = 4, fig.height = h)
#+ cellcylce_pre, include=TRUE
gp_pre

opts_chunk$set(fig.width= w, fig.height = h)
#+ cellcylce_post, include=TRUE
gp_post


## ######################################################################### ##
## ####################### (v) inspect PCA ################################# ##
## ######################################################################### ##

gp1 <- DimPlot(object = s.full, reduction = "pca", pt.size = .1, 
              group.by = opt$split_var)
gp1_2 <- DimPlot(object = s.full, reduction = "pca", pt.size = .1, dims = c(1,3),
               group.by = opt$split_var)
gp1_3 <- DimPlot(object = s.full, reduction = "pca", pt.size = .1, dims = c(2,3),
                 group.by = opt$split_var)
gp2 <- ElbowPlot(s.full, reduction = "pca")
gp3 <- VizDimLoadings(s.full, dims = 1:2, reduction = "pca")


#' ## Inspection of PCA
#' The plots below show the first principle components colored by `r opt$split_var`,
#'  the elbow plot and the genes with highest loadings.
opts_chunk$set(fig.width = 8, fig.height = 4)
#+ pca1, include=TRUE
grid.arrange(gp1, gp1_2, ncol=2)
opts_chunk$set(fig.width = 8, fig.height = 4)
#+ pca2, include=TRUE
grid.arrange(gp1_3, gp2, ncol = 2)

opts_chunk$set(fig.width = 8, fig.height = 4)
#+ pca3, include=TRUE
gp3


pca_coord <- as.data.frame(Embeddings(s.full, reduction = "pca")[,1:12])
pca_coord$sample <- s.full[[opt$split_var]][,opt$split_var]

melted <- melt(pca_coord, id.vars = "sample")

gp <- ggplot(melted, aes(x=value, color = sample)) + geom_density()
gp <- gp + facet_wrap(~variable, ncol = 4) + theme_bw()
if (length(unique(melted$sample)) > 10 ){
  ## do not include legend if too many conditions/samples
  gp <- gp + theme(legend.position = "none")
} else {
  #gp <- gp + guides(col = guide_legend(ncol = 2))
  gp <- gp + theme(legend.text = element_text(size = 8))  
}


#' ## Inspection of difference in PCA across samples/conditions
#' The plots below show density plots of the distribution of values for each PC
#' in the different samples. If more than 10 samples are included, the legend is not 
#' shown for space reasons.
opts_chunk$set(fig.width = 15, fig.height = 13)
#+ pca_density, include=TRUE
gp

flog.info("Finished PCA plots")
flog.info("Completed")
