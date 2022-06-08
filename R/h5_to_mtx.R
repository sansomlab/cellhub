stopifnot(
  require(optparse),
  require(futile.logger),
  require(BiocParallel),
  require(Matrix),
  require(Seurat),
  require(SeuratDisk),
  require(dplyr),
  require(DropletUtils) 
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--h5_file"), default="sc.h5",
              help="*.h5 file (e.g. cellbender output)"),
  make_option(c("--features"), default="features.tsv.gz",
              help="Single column file features.tsv.gz with the gene names
              of .h5 file"),
  make_option(c("--sample_id"), default="sample1",
              help="String with sample name to add to barcode IDs."),
  make_option(c("--numcores"), default=1,
              help="Number of cores."),
  make_option(c("--log_filename"), default="h5_to_mtx.log"),
  make_option(c("--output_dir"), default="./mtx")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger -----------------------------------------------------------------------
# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

# Number of cores to use -------------------------------------------------------

multicoreParam <- MulticoreParam(workers = opt$numcores)

# Read in data -----------------------------------------------------------------

flog.info("Reading h5 object...")

if(file.exists(opt$h5_file)) {
  h5 <- Read10X_h5(filename = opt$h5_file, use.names = FALSE)
  flog.info("h5 object: ", str(h5), capture = TRUE)
} else {
  stop("Provided h5 file does not exits.")
}

# create Seurat object ---------------------------------------------------------
flog.info("Creating seurat object...")

obj <- CreateSeuratObject(counts = h5)
flog.info("seurat object: ", obj, capture = TRUE)
flog.info("genes: ", head(row.names(obj)), capture = TRUE)
flog.info("barcodes: ", head(colnames(obj)), capture = TRUE)

if(file.exists(opt$features)) {
  
  feat <- read.delim(opt$features, header=FALSE)
  flog.info("feature file: ", head(feat), capture = TRUE)
  
  if(all(rownames(obj) == feat$V1)) {
    
    barcode_ids <- paste0(gsub("-.+", "", colnames(obj)),
                          "-", opt$sample_id)
    
    flog.info("barcode names: ", head(barcode_ids), capture = TRUE)
    
    write10xCounts(path = opt$output_dir,
                   obj@assays$RNA@counts,
                   gene.id=row.names(obj),
                   gene.symbol=feat$V2,
                   barcodes=barcode_ids,
                   version="3",
                   overwrite = TRUE)
    
    flog.info("mtx matrices written in: ", opt$output_dir, capture = TRUE)
    
  } else {
    stop("Provided feature file does not contain h5 features in first column.")
  }
} else {
  stop("Provided feature file does not exits.")
}