# Libraries ------
stopifnot(require(optparse),
          require(futile.logger),
          require(dplyr),
          require(DropletUtils),
          require(tibble),
          require(scater),
          require(Matrix),
	  require(BiocParallel)
#	  require(cellqc)
)

# modify from https://stackoverflow.com/questions/27418461
find_valley <- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] < x[i-1]) & (x[i] < x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}

# Global options ------

options(stringsAsFactors = F)


# Script arguments ------

option_list <- list(
  make_option(c("--cellranger_dir"), default=".",
              help="Folder with filtered cellranger output. Must include barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--library_id"), default=NULL,
              help="library or channel id"),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to run scater's perCellQCMetrics function"),
  make_option(c("--log_filename"), default="qcmetrics.log"),
  make_option(c("--outfile"), default="qcmetrics.tsv.gz")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger ------

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

# Number of cores to use perCellQC function ------

multicoreParam <- MulticoreParam(workers = opt$numcores)


# Read in data ------
#####################

flog.info("Reading data and making SingleCellExperiment object...")
s <- read10xCounts(file.path(opt$cellranger_dir))
flog.info("Number of cells in input: %s", format(ncol(s), big.mark=","))


# Compute basic QC metrics ------
#################################

flog.info("Calculating basic QC metrics...")

# Number of genes
ngenes <- colSums(counts(s) > 0)

# Total UMI counts
total_UMI <- colSums(counts(s))

# Create dataframe
cell_qc <- tibble(barcode_id = colData(s)$Barcode,
                  library_id = opt$library_id,
                  ngenes = ngenes,
                  total_UMI = total_UMI,
                  log_total_UMI = log10(total_UMI))

adj <- 1
limits <- NULL
while(length(limits) != 2) {
  den <- density(filter(cell_qc, log_total_UMI > 0)$log_total_UMI, adjust = adj)
  adj <- adj + 0.05
  d <- den$y
  log_depth <- den$x
  limits <- log_depth[find_valley(d)]
  print(limits)
}

if(length(limits) == 2) {
  
  if(limits[2] > 2 & limits[2] < 3) {
    
    # background limits
    blim_min <- limits[1] + 0.1
    blim_max <- limits[2] - 0.1
    # drops limit
    dlim_min <- limits[2] + 0.1
    cell_qc %>%
      mutate(group = ifelse(log_total_UMI < blim_min, "debri", 
                            ifelse(log_total_UMI > blim_min & log_total_UMI < blim_max, "background",
                                   ifelse(log_total_UMI > dlim_min, "cell", "intermediate")))) -> cell_qc
    
    flog.info(table(cell_qc$group))
              
  } else if(limits[1] > 2 & limits[1] < 3) {
    
    # background limits
    blim_max <- limits[1] - 0.1
    # drops limit
    dlim_min <- limits[1] + 0.1
    cell_qc %>%
      mutate(group = ifelse(log_total_UMI < blim_max, "background", 
                            ifelse(log_total_UMI > dlim_min, "cell", "intermediate"))) -> cell_qc
    
    flog.info(table(cell_qc$group))
              
  }

} else {
  
  cell_qc %>%
    mutate(group = "undetermined") -> cell_qc
  
}

# Write table
flog.info("Writing output table")
write.table(cell_qc, file = gzfile(opt$outfile),
            quote=FALSE, row.names = FALSE, sep="\t")

flog.info("Completed")