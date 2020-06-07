## Title ----
##
## Run harmony
##
## Description ----
## Run integration using harmony
##
# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(Matrix))
stopifnot(require(harmony))
stopifnot(require(knitr))
stopifnot(require(future))
stopifnot(require(futile.logger))

# Parameters -------------------------------------------------------------------
# The script expects the following parameters:

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

# these options are read from the yml
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
  # cell cycle already added to this as part of the pipeline.py
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

## add cell cycle to regression variables
if ( c("all") %in% opt$regress_latentvars ) {
  flog.info("Cell cycle correction for S and G2M scores will be applied")
  vreg <- c(vreg[!grepl("all", vreg)], "S.Score", "G2M.Score")
} else if ( c("difference") %in% opt$regress_latentvars ){
  flog.info("Cell cycle correction for the difference between G2M and S phase scores will be applied")
  vreg <- c(vreg[!grepl("difference", vreg)], "CC.Difference")
} else {
  flog.info("Data will be scaled without correcting for cell cycle")
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


## ######################################################################### ##
## ####################### (iii) integrate data ############################ ##
## ######################################################################### ##

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

flog.info("Default Assay after integration: %s", DefaultAssay(replicates.integrated))
assay <- DefaultAssay(replicates.integrated)


## write out the integrated coordinates (non-UMAP)
int_out <- data.frame(Embeddings(object = replicates.integrated,
                                 reduction = integrated_reduction))
int_out$barcode <- rownames(int_out)

write.table(int_out, gzfile(file.path(opt$outdir, paste0(opt$tool,".tsv.gz"))),
            sep="\t", quote=FALSE, row.names=FALSE)

# write out the stdev slot for plotting
if (opt$tool %in% c("rawdata")){
  flog.info("No integration performed, so no stdev to plot")
} else {
  if (length(replicates.integrated@reductions$harmony@stdev) != 0) {
    df = data.frame(comp = c(1:length(replicates.integrated@reductions$harmony@stdev)), 
                  stdev = replicates.integrated@reductions$harmony@stdev)
    write.table(df, gzfile(file.path(opt$outdir, paste0(opt$tool,"_stdev.tsv.gz"))),
              sep="\t", quote=FALSE, row.names=FALSE)
  }
}

flog.info("Finished integration and saved reduced dimension coordinates ...")

## ######################################################################### ##
## ##################### (iv) Run UMAP ##################################### ##
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
plot_data <- merge(umap, metadata, by="barcode")

write.table(plot_data, gzfile(file.path(opt$outdir, "umap.tsv.gz")),
            sep="\t", quote=FALSE, row.names=FALSE)

flog.info("Finished UMAP and saving UMAP coordinates ...")

## ######################################################################### ##
## ################### (v) Save integrated data ############################ ##
## ######################################################################### ##

# Add s.full_misc
replicates.integrated@misc <- s.full_misc

flog.info("Saving object ...")
flog.info("Outfile: %s", file.path(opt$outdir, "integrated.rds"))
saveRDS(replicates.integrated, file = file.path(opt$outdir, "integrated.rds"))

flog.info("Completed")



