#' ---
#' title: "Run log-normalization and assessment of highly variable genes"
#' output:
#'  html_document:
#'   self_contained: false
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
stopifnot(require(tenxutils))
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
  # normalized and scaled already
  "seurat_obj" = "",

  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",

  # (Harmony) Whether to perform normalisation per batch or on merged
  # object. By default is set to 'merged' according to current tutorial.
  # merged | perbatch
  "merge_normalisation" = "merged",
  
  # Model used to regress out latent variables (log-normalisation)
  "regress_modeluse" = "linear",

  # Variables to regress, example: percent.mito
  "regress_latentvars" = NULL,

  # Whether to regress cell cycle, options: none, difference, all
  "regress_cellcycle" = "none",

  # A vector of Ensembl gene ids associated with S phases.
  # See Seurat::CellCycleScoring(s.genes=...)
  "sgenes" = NULL,

  # A vector of Ensembl gene ids associated with G2M phase.
  # See Seurat::CellCycleScoring(g2m.genes=...)
  "g2mgenes" = NULL,

  # Number of cores to use
  "numcores" = 1,
  
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

plan("multiprocess",
     workers = opt$numcores)

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

flog.info("Reading RDS file...  ")
s.full <- readRDS(opt$seurat_obj)
flog.info("Default assay of Seurat object: %s", DefaultAssay(s.full))

# Make a dataframe of s@misc to add is to the final seurat object
s.full_misc <- s.full@misc

# # reset to RNA as the SCT transform was only an initial one!
# DefaultAssay(s.full) <- "RNA"
# flog.info("Reset default assay of Seurat object: %s", DefaultAssay(s.full))


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
  sgenes <- read.table(opt$sgenes, header=F, as.is=T)$V1
  sgenes <- s.full@misc$seurat_id[s.full@misc$gene_id %in% sgenes]

  g2mgenes <- read.table(opt$g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s.full@misc$seurat_id[s.full@misc$gene_id %in% g2mgenes]

  s.full <- RunPCA(object = s.full, features = c(sgenes, g2mgenes), verbose = FALSE,
                   reduction.name = "pcacellcyclePre", reduction.key="pcacellcyclePre_")
  gp_pre <- DimPlot(object = s.full, group.by="Phase", reduction = "pcacellcyclePre")
}

if (opt$merge_normalisation %in% c("merged")){
  flog.info("Running on merged object, so Seurat object will not be split ...  ")
  txt_merged <- "Normalization was performed on the merged object."
} else {
  txt_merged <- "Normalization was performed per batch."
  s.list <- SplitObject(s.full,
                      split.by = opt$split_var)
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

flog.info("Final variables for regression are %s", paste0(vreg, collapse = ", "))

## ######################################################################### ##
## ################ (ii) Normalize data and find hvg ####################### ##
## ######################################################################### ##

#' ## Data normalization
#'
#' The data was normalized using log-normalization.
#' `r txt_merged`
#' The following variables were regressed: `r paste(vreg, collapse=",")` .

flog.info("Performing log-normalization ...")

if (opt$merge_normalisation %in% c("merged")) {
  ## Perform log-normalization of the RNA assay
  flog.info("Merged log-normalization ...")
  s.full <- NormalizeData(object=s.full,
                          normalization.method="LogNormalize",
                          scale.factor=10E3, verbose = FALSE)
  
  s.full <- FindVariableFeatures(s.full,
                                 selection.method = "vst",
                                 nfeatures = opt$ngenes, 
                                 verbose = TRUE) 
} else if (opt$merge_normalisation %in%  c("perbatch")) {
  flog.info("Per-batch log-normalization ...")
  for (i in 1:length(s.list)) {
    ## Perform log-normalization of the RNA assay
    s.list[[i]] <- NormalizeData(object=s.list[[i]],
                                 normalization.method="LogNormalize",
                                 scale.factor=10E3, verbose = FALSE)
    
    s.list[[i]] <- FindVariableFeatures(s.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = opt$ngenes, 
                                        verbose = FALSE)

  }
}

flog.info("Done log-normalization for samples ...")


## ######################################################################### ##
## ###################### (iii) assess hv genes ############################ ##
## ######################################################################### ##

flog.info("Assess highly variable genes ...")
## Assessment of highly-variable genes
gglist <- list()
hvg_info <- data.frame()

xvar = "mean"
yvar = "variance.standardized"


if (opt$merge_normalisation %in% c("merged")) {
  max_iter = 1
  loop = FALSE
} else { 
  max_iter = length(s.list)
  loop = TRUE
}
  
for (i in 1:max_iter) {
  if (loop){
    n <- names(s.list)[i]
    hvg <- HVFInfo(object = s.list[[i]] )
    hvg$var.gene = FALSE
    hvg$var.gene[rownames(hvg) %in% VariableFeatures(object = s.list[[i]])] <- TRUE
    hvg$sample = n
    hvg_info <- rbind(hvg_info, hvg) 
  } else {
    n <- "hvg_merged"
    hvg <- HVFInfo(object = s.full)
    hvg$var.gene = FALSE
    hvg$var.gene[rownames(hvg) %in% VariableFeatures(object = s.full)] <- TRUE
    hvg$sample = n
    hvg_info = hvg
  }
  hvg$dispersion = 0.1
  
  melted <- melt(hvg[, c(xvar, yvar, "var.gene")],
                 id.vars=c("var.gene", xvar))
  
  gp <- ggplot(melted, aes_string(xvar, "value", color="var.gene"))
  gp <- gp + scale_color_manual(values=c("black","red"))
  gp <- gp + geom_point(alpha = 1, size=0.5)
  gp <- gp + theme_bw()
  gp <- gp + ylab(yvar) + ggtitle(n)
  gp <- gp + theme(plot.title = element_text(hjust = 0.5))
  
  if (loop){
    gglist[[n]] <- gp
  } else {
    break
  }
}

if(loop) {
  nr <- ceiling(length(s.list)/3)
  gp <- plot_grid(plotlist = gglist, 
                  nrow=nr)
}

w <- min(3, max_iter)*4 + 1
h <- max(1, max_iter/3)*4

save_ggplots(file.path(opt$outdir, "var_genes_plot"),
             gp=gp,
             width=w+1,
             height=h)

#' ## Assessment of highly-variable genes
#' Number of variable genes selected per sample: `r opt$ngenes` genes.
opts_chunk$set(fig.width= w, fig.height = h)
#+ hvg_genes, include=TRUE
gp

flog.info("Finished plots for hv genes ...")

# add gene_id to table
hvg_info$gene_name <- rownames(hvg_info)
hvg_info$gene_id <- s.full_misc$gene_id[match(hvg_info$gene_name, s.full_misc$seurat_id)]

flog.info("Keeping only genes that are variable in any sample ...")
hvg_info = hvg_info[hvg_info$var.gene == TRUE,]

# add info from genelists that were specified
if (!is.null(opt$vargenes_dir)){
  flog.info("Check hv genes for given genelists ...")
  names_genesets = c()
  for (f in list.files(opt$vargenes_dir, pattern = "*.tsv*")){
    if (endsWith(f, ".gz")){
      list_name = substr(f, 1, nchar(f)-nchar(".tsv.gz"))
      } else {
      list_name = substr(f, 1, nchar(f)-nchar(".tsv"))
      }
    names_genesets = c(names_genesets, list_name)
    genelist = read.table(file.path(opt$vargenes_dir, f), header = TRUE)
    hvg_info[, list_name] <- FALSE
    hvg_info[hvg_info$gene_id %in% genelist$gene_id, list_name] <- TRUE
  }
  flog.info("Make plots of fraction of hvg within each given geneset.")
  df_plot = hvg_info[, c("sample", names_genesets)] %>% pivot_longer(all_of(names_genesets)) %>%
            group_by(sample, name) %>% summarise(n_set = n(),
                                                 fraction_geneset = length(which(value == TRUE))/n_set) 
  nsamples = length(unique(df_plot$sample))
  if(opt$merge_normalisation %in% c("perbatch")){ 
    df_plot$sample = factor(as.character(df_plot$sample), levels = names(s.list))
    }
  gp <- ggplot(df_plot, aes(x=name, y=fraction_geneset)) + geom_bar(stat='identity', fill = "grey")
  gp <- gp + facet_wrap(~sample, ncol = min(nsamples, 3)) + ylab("Fraction of geneset \n within highly variable genes")
  gp <- gp + xlab("Geneset") + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                                                  strip.background = element_blank())
  #gp <- gp + coord_cartesian(ylim = c(0,1))
  
  w <- min(3, max_iter)*4
  h <- max(1, max_iter/3)*4 + 2
  
  save_ggplots(file.path(opt$outdir, "hvg_fraction_genelists"),
               gp=gp,
               width=w,
               height=h)
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



## ######################################################################### ##
## ####################### (iii) scale data ################################ ##
## ######################################################################### ##


flog.info("Run data scaling...")

if (opt$merge_normalisation %in% c("merged")) {
  ## Scaling of the RNA assay data
  #all.genes <- rownames(s.full)
  s.full <- ScaleData(object=s.full,
                      #features = all.genes,
                      vars.to.regress=vreg,
                      model.use=opt$regress_modeluse)
  flog.info("Running PCA ...")
  s.full <- RunPCA(s.full, npcs = 50, verbose = FALSE)
} else {
  for (i in 1:length(s.list)) {
  ## Scaling of the RNA assay data
  #all.genes <- rownames(s.list[[i]])
  s.list[[i]] <- ScaleData(object=s.list[[i]],
                           #features = all.genes,
                           vars.to.regress=vreg,
                           model.use=opt$regress_modeluse)
  }
  features <- SelectIntegrationFeatures(object.list = s.list, 
                                        nfeatures = opt$ngenes)
  flog.info("Merge Seurat objects...")
  s.full <- s.list[[1]]
  for (i in 2:length(s.list)){
    s.full <- merge(s.full, s.list[[i]], merge.data = TRUE)
  }
  flog.info("Re-scaling merged data ...")
  # ## Scaling of the RNA assay data
  s.full <- ScaleData(object=s.full,
                      features = features,
                      vars.to.regress=vreg,
                      model.use=opt$regress_modeluse, verbose = FALSE)
  s.full <- RunPCA(s.full, npcs = 50, features = features,
                   verbose = FALSE)
}


if (exists("features")){
  flog.info("Add info whether gene was used for integration to hv table ...")
  hvg_info$use_integration = FALSE
  hvg_info$use_integration[hvg_info$gene_name %in% features] = TRUE
} else {
  flog.info("All genes are used for integration.")
  hvg_info$use_integration = TRUE
}

write.csv(hvg_info, gzfile(file.path(opt$outdir, "hv_genes_info.csv.gz")), 
          quote = F, row.names = F)
flog.info("Finished writing file with hv genes ...")


#### check cell cycle regression

if ( identical(opt$regress_cellcycle, "none") ) {
  flog.info("No cell cycle plots produced post correction...")
  gp_post <- ggplot() + theme_void()
  h = 0.1
  w = 4
} else { 
  if (opt$merge_normalisation %in% c("merged")){
    s.full <- RunPCA(object = s.full, features = c(sgenes, g2mgenes), 
                     verbose = FALSE, reduction.name = "pcacellcyclePost", 
                     reduction.key="pcacellcyclePost_")
    gp_post <- DimPlot(object = s.full, group.by="Phase",
                       reduction = "pcacellcyclePost")
    w=4
  } else {
    gglist <- list()
    for (i in 1:length(s.list)) {
      s.list[[i]] <- RunPCA(object = s.list[[i]], features = c(sgenes, g2mgenes),
                            verbose = FALSE, reduction.name = "pcacellcyclePost",
                            reduction.key="pcacellcyclePost_")
      gglist[[i]] <- DimPlot(object = s.list[[i]], group.by="Phase",
                             reduction = "pcacellcyclePost")
    }
    nr <- ceiling(length(s.list)/3)
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


write.table(Embeddings(s.full, reduction = "pca"), 
            file = gzfile(file.path(opt$outdir, "pca_comp.tsv.gz")))

write.table(Loadings(s.full, reduction = "pca"), 
            file = gzfile(file.path(opt$outdir, "pca_loadings.tsv.gz")))

flog.info("Finished PCA plots and writing out pca components")


# Add s.full_misc
if (opt$merge_normalisation %in% c("perbatch"))
s.full@misc <- s.full_misc

flog.info("Saving object ...")
flog.info("Outfile: %s", file.path(opt$outdir, "pre_integrated.rds"))
saveRDS(s.full, file = file.path(opt$outdir, "pre_integrated.rds"))

flog.info("Completed")
