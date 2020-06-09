## Title ----
##
## Create Seurat object
##
## Description ----
## Create a Seurat object
##
## Default parameters ----
## example yml: /gfs/devel/kjansen/dropflow/Rmd/emptydrops.yml

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(Matrix))
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
  # The path to the directory containing matrices and associates
  # barcodes, features and metadata
  "matrixdir" = "",

  # The path to the Seurat object (output)
  "seurat_obj" = "",

  # The path to the output directory
  "outdir" = "",

  # The name of the project
  "project" = "",

  # Latent variables to regress out, ONLY used for cell cycle scoring
  # See Seurat::ScaleData(vars.to.regress=..., model.use=opt$modeluse)
  "latentvars" = "none",

  # Model used to regress out latent variables
  "modeluse" = "linear",

  # A vector of Ensembl gene ids associated with S phases.
  # See Seurat::CellCycleScoring(s.genes=...)
  "sgenes" = "none",

  # A vector of Ensembl gene ids associated with G2M phase.
  # See Seurat::CellCycleScoring(g2m.genes=...)
  "g2mgenes" = "none",

  # Number of cores to use
  "numcores" = 1
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

flog.info("Running with options: ", opt, capture = TRUE)
flog.info("\n")

plan("multiprocess",
     workers = opt$numcores)

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

flog.info("Importing matrix from: %s", opt$matrixdir)
data <- Read10X(opt$matrixdir)
flog.info("Finished.")

## Seurat discards the Ensembl IDs and makes it's own identifiers (!)
## From https://github.com/satijalab/seurat/blob/master/R/preprocessing.R:
##
## rownames(data) <- make.unique(
##    names=as.character(
##        x=sapply(
##            X=gene.names,
##            FUN=ExtractField,
##            field=2,
##            delim="\\t"
##
## We need to track the seurat id -> Ensembl id mapping, e.g. for downstream GO analysis
## (and to guarentee reproducibility, e.g. between annotations versions)

inFile <- file.path(opt$matrixdir, "features.tsv.gz")
flog.info("Importing gene information from: %s", inFile)
genes <- read.table(gzfile(inFile), as.is=TRUE)
genes <- genes[,c(1,2)]

colnames(genes) <- c("gene_id", "gene_name")
genes$seurat_id <- as.character(make.unique(genes$gene_name))
rownames(genes) <- genes$seurat_id

write.table(
  genes, gzfile(file.path(opt$outdir, "annotation.txt.gz")),
  sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)

## ######################################################################### ##
## ################## (ii) Create Seurat object ############################ ##
## ######################################################################### ##

flog.info("Creating Seurat object ... ")
s <- CreateSeuratObject(counts=data,
                        project=opt$project
)
flog.info("Finished.")

s@misc <- genes

## ######################################################################### ##
## ################## (ii) Add metadata #################################### ##
## ######################################################################### ##

metaFile <- file.path(opt$matrixdir, "metadata.tsv.gz")

## Read in the metadata
if (!file.exists(metaFile)) {
  stop("No metadata file given")
}

flog.info("Importing metadata ... ")
metadata <- read.table(gzfile(metaFile),
                       sep="\t", header=TRUE, as.is=TRUE)

# ensure that the barcode column is present in the metadata
if (!"barcode" %in% colnames(metadata)) {
  stop('Mandatory "barcode" column missing from the metadata')
}

rownames(metadata) <- metadata$barcode
metadata$barcode <- NULL

metadata <- metadata[colnames(x = s), ]

flog.info("Adding meta data ... ")
for(meta_col in colnames(metadata))
{
  s[[meta_col]] <- metadata[[meta_col]]
}
flog.info("Finished adding meta data")


## ######################################################################### ##
## ######### (vi) Optional calculation of cell cycle scores ################ ##
## ######################################################################### ##

## Cell cycle correction requires initial normalisation and scaling, but no
## cell cycle correction is applied at this stage.

if(identical(opt$latentvars, "none")) {
  latent.vars = NULL
  flog.info("No latent vars specified")
} else {
  if(grepl(",", opt$latentvars)){
    latent.vars <- unlist(strsplit(opt$latentvars, split=","))
  } else {
    latent.vars <- opt$latentvars
  }
  flog.info("Latent vars for initial normalisation: %s", latent.vars)
}


flog.info("Performing initial log-normalization without cell cycle")

## Perform log-normalization of the RNA assay
s <- NormalizeData(object=s,
                   normalization.method="LogNormalize",
                   scale.factor=10E3)

## Initial scaling of the RNA assay data
all.genes <- rownames(s)
s <- ScaleData(object=s,
               features = all.genes,
               vars.to.regress=latent.vars,
               model.use=opt$modeluse)

flog.info("Done log-normalization")

## If cell cycle genes are given, make a PCA of the cells based
## on expression of cell cycle genes

if (!(is.null(opt$sgenes) | opt$sgenes=="none")
    & !(is.null(opt$g2mgenes) | opt$g2mgenes=="none")) {

  flog.info("Cell cycle genes are present...")

  # get the genes representing the cell cycle phases
  sgenes_ensembl <- read.table(opt$sgenes, header=F, as.is=T)$V1
  sgenes <- s@misc$seurat_id[s@misc$gene_id %in% sgenes_ensembl]

  g2mgenes_ensembl <- read.table(opt$g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s@misc$seurat_id[s@misc$gene_id %in% g2mgenes_ensembl]

  # score the cell cycle phases
  s <- CellCycleScoring(object=s,
                        s.features=sgenes,
                        g2m.features=g2mgenes,
                        set.ident=TRUE)

  s$CC.Difference <- s$S.Score - s$G2M.Score

  flog.info("Finished cell cycle scoring and added metadata column...")

}

## ######################################################################### ##
## ################ (vii) Save Seurat object ############################### ##
## ######################################################################### ##

flog.info("seurat object final default assay: %s", DefaultAssay(s))

# Save the R object
saveRDS(s, file=opt$seurat_obj)

flog.info("Completed")
