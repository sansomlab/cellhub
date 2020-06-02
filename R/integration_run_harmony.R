#' ---
#' title: "Run harmony integration"
#' output:
#'  html_document:
#'   self_contained: false
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
stopifnot(require(Seurat))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(Matrix))
stopifnot(require(future))
stopifnot(require(SeuratWrappers))
stopifnot(require(harmony))
stopifnot(require(tenxutils))
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
  
  # Path to the file with highly variable features that should be used for
  # the PCA + harmony integration
  "hv_genes" = NULL,
  
  
  # Variables to regress, example: percent.mito
  "regress_latentvars" = NULL,
  
  # Model used to regress out latent variables (log-normalisation)
  "regress_modeluse" = "linear",
  
  # Harmony or rawdata
  "tool" = "harmony",
  
  # Variable for harmony, 30 used for rawdata
  "nPCs" = 30,  
  
  # (Harmony) Width of soft kmeans clusters, default = 0.1
  "sigma" = 0.1,
  
  # (Harmony) Diversity clustering penalty parameter, default = 2
  "theta" = 2,
  
  # (Harmony) Ridge regression penalty parameter, default = 1
  "lambda" = 1,
  
  # Number of cores to use
  "numcores" = 1
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

## ######################################################################### ##
## ################ (ii) re-run PCA on selected hv genes ################### ##
## ######################################################################### ##

if (endsWith(opt$hv_genes, ".gz")) {
  hvg_info <- read.csv(gzfile(opt$hv_genes), header = TRUE)
} else {
  hvg_info <- read.csv(opt$hv_genes, header = TRUE)
}

features <- hvg_info$gene_name[hvg_info$use_integration == TRUE]

if(! all(features %in% rownames(s.full))) {
  stop("Some features provided as highly variable genes are absent from the
       seurat object")
}

assay_use = "RNA"

flog.info("Rescale the data using the input hv genes...")
# rescale here in case the hv genes have changed.
s.full <- ScaleData(object=s.full,
                    features = features,
                    vars.to.regress=vreg,
                    model.use=opt$regress_modeluse)
# reset the variable features instead of inputting them directly into PCA
# (otherwise error regarding subscript out of bounds)
VariableFeatures(s.full) <- features
s.full <- RunPCA(s.full, npcs = opt$nPCs, #features = features,
                 verbose = FALSE)

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
## ####################### (iii) integrate data ############################ ##
## ######################################################################### ##

#' ## Data integration using: `r opt$tool`

# unclear if harmony dims.use argument is working, therefore components defined
# above in RunPCA statement and then full 'pca' slot passed to harmony.
# issue here: https://github.com/immunogenomics/harmony/issues/82

if (opt$tool == "harmony"){
  replicates.integrated <- RunHarmony(s.full,
                                      group.by.vars = opt$split_var,
                                      verbose = FALSE,
                                      lambda = opt$lambda,
                                      theta = opt$theta, sigma = opt$sigma,
                                      assay.use = assay_use)
  integrated_reduction = "harmony"
  nPCs <- min(dim(replicates.integrated@reductions$harmony@cell.embeddings)[2],50)
  flog.info("Done integration using harmony ...") 
} else if (opt$tool == "rawdata") {
  flog.info("Skipped integration ...")
  replicates.integrated <- s.full
  integrated_reduction = "pca"
}

if (opt$tool %in% c("harmony")){
  features_used = paste0("The integration was performed using ", 
                         opt$nPCs," PCs for harmony")
} else if (opt$tool == "rawdata") {
  features_used = "No integration was performed"
}

#' `r features_used`.
#'

flog.info("Default Assay after integration: %s", DefaultAssay(replicates.integrated))
assay <- DefaultAssay(replicates.integrated)

## ######################################################################### ##
## ############### (iv) Downstream analyses and plots ###################### ##
## ######################################################################### ##


# Run UMAP
flog.info("Running UMAP ...")
n_comp = 30
  
replicates.integrated <- RunUMAP(replicates.integrated,
                                 reduction = integrated_reduction,
                                 dims = 1:n_comp,
                                 assay = assay)
replicates.integrated@meta.data$barcode = rownames(replicates.integrated[[]])

## extract the UMAP coordinates from the seurat object
umap <- data.frame(Embeddings(object = replicates.integrated,
                              reduction = "umap"))
umap$barcode <- rownames(umap)
metadata = replicates.integrated[[]]
metadata$barcode = rownames(metadata)
plot_data <- left_join(umap, metadata, by="barcode")

##### UMAP PLOTS
gp <- ggplot(plot_data, aes_string(x="UMAP_1", y="UMAP_2", col=opt$split_var))
gp <- gp +  geom_point(alpha=0.5, size=2) + theme_bw()
gp <- gp + guides(col = guide_legend(override.aes = list(alpha=1)))

save_ggplots(file.path(opt$outdir, "integrated_umap_nofacet"),
             gp)

#' UMAP (`r n_comp` components) colored by variable that was used for integration.
#+ umap, include=TRUE, fig.height=4, fig.width=6
print(gp)

gp <- gp + facet_wrap(as.formula(paste("~", opt$split_var)), ncol = 3)


w <- min(3, length(unique(s.full[[]][, opt$split_var])))*4
h <- max(1, length(unique(s.full[[]][, opt$split_var]))/3)*4

save_ggplots(file.path(opt$outdir, "integrated_umap_facet_splitvar"),
             gp, width = w, height = h)

opts_chunk$set(fig.width= w, fig.height = h)
#' UMAP (`r n_comp` components) split by the variable that was used for integration.
#+ umap_facet, include=TRUE
print(gp)
#+ include=FALSE

write.table(plot_data, file.path(opt$outdir, "umap.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

if (opt$tool %in% c("harmony")) {
  if (length(replicates.integrated@reductions$harmony@stdev) == 0) {
    flog.info("cannot make elbow plot from harmony dimensions because stdev is absent.")
    h=1
    w=6
    elbow_plot <- ggplot() + theme_void()
    txt <- "No standard deviation was present, therefore no elbow plots."
  } else {
    h=4
    w=4
    elbow_plot <- ElbowPlot(replicates.integrated, ndims=nPCs, reduction = "harmony")
    txt <- "Elbow plot for harmony stdev shown below."
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


## ######################################################################### ##
## ################### (v) Save integrated data ############################ ##
## ######################################################################### ##

# Add s.full_misc
replicates.integrated@misc <- s.full_misc

flog.info("Saving object ...")
flog.info("Outfile: %s", file.path(opt$outdir, "integrated.rds"))
saveRDS(replicates.integrated, file = file.path(opt$outdir, "integrated.rds"))

## write out the integrated coordinates (non-UMAP)
int_out <- data.frame(Embeddings(object = replicates.integrated,
                                 reduction = integrated_reduction))
int_out$barcode <- rownames(int_out)
metadata = replicates.integrated[[]]
metadata$barcode = rownames(metadata)
plot_data <- left_join(int_out, metadata, by="barcode")

write.table(plot_data, gzfile(file.path(opt$outdir, paste0(opt$tool,".tsv.gz"))),
            sep="\t", quote=FALSE, row.names=FALSE)

flog.info("Completed")



