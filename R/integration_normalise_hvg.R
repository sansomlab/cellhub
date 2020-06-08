## Title ----
##
## Normalisation and highly-variable genes
##
## Description ----
## Run normalisation and determine highly-variable genes


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

  # Path to the seurat begin.rds object
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
  "regress_latentvars" = "none",

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
if ( ! identical(opt$regress_latentvars, "none")){
  if(grepl(",", opt$regress_latentvars)){
    vreg <- unlist(strsplit(opt$regress_latentvars, split=","))
  } else {
    vreg <- opt$regress_latentvars
  }
} else {
  vreg = NULL
  flog.info("No variable regression...")
}


## run PCA pre normalisation if cell cycle is included.

if ( identical(opt$regress_cellcycle, "none") ) {
  flog.info("No plots produced")
} else {
  sgenes <- read.table(opt$sgenes, header=F, as.is=T)$V1
  sgenes <- s.full@misc$seurat_id[s.full@misc$gene_id %in% sgenes]

  g2mgenes <- read.table(opt$g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s.full@misc$seurat_id[s.full@misc$gene_id %in% g2mgenes]

  s.full <- RunPCA(object = s.full, features = c(sgenes, g2mgenes), verbose = FALSE,
                   reduction.name = "pcacellcyclePre", reduction.key="pcacellcyclePre_")
  flog.info("Create and save the plots here as the slot after splitting of the Seurat
            object...")
  gp_pre <- DimPlot(object = s.full, group.by="Phase", reduction = "pcacellcyclePre")
  saveRDS(gp_pre, file.path(opt$outdir, "gp_cellcyclePre"))
}

if (opt$merge_normalisation %in% c("merged")){
  flog.info("Running on merged object, so Seurat object will not be split...")
} else {
  flog.info("Normalization was performed per batch.")
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


## ######################################################################### ##
## ################ (ii) Normalize data and find hvg ####################### ##
## ######################################################################### ##

flog.info("Final variables for regression are %s", paste0(vreg, collapse = ", "))
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
    gglist <- gp
    break
  }
}

saveRDS(gglist, file.path(opt$outdir, "gglist_hvg.rds"))


flog.info("Finished assessment of hv genes, saved plots ...")

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
} else {
  flog.info("No genelists provided")
}

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
  ## Scaling of the RNA assay data
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
} else { 
  if (opt$merge_normalisation %in% c("merged")){
    s.full <- RunPCA(object = s.full, features = c(sgenes, g2mgenes), 
                     verbose = FALSE, reduction.name = "pcacellcyclePost", 
                     reduction.key="pcacellcyclePost_")
  } else {
    flog.info("Post cell cycle correction PER batch is saved as information is lost
              in Seurat object, that was re-merged...")
    gglist <- list()
    for (i in 1:length(s.list)) {
      s.list[[i]] <- RunPCA(object = s.list[[i]], features = c(sgenes, g2mgenes),
                            verbose = FALSE, reduction.name = "pcacellcyclePost",
                            reduction.key="pcacellcyclePost_")
      gglist[[i]] <- DimPlot(object = s.list[[i]], group.by="Phase",
                             reduction = "pcacellcyclePost")
    }
    saveRDS(gglist, file.path(opt$outdir, "gglist_cellcyclePost.rds"))
  }
}

## ######################################################################### ##
## ####################### (v) write out PCA info ########################## ##
## ######################################################################### ##


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
