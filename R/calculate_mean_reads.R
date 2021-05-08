## Title ----
##
## Mean UMI per cell
##
## Description ----
## Calculate the mean total UMI per cell for all barcodes

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(DropletUtils))
stopifnot(require(futile.logger))
stopifnot(require(knitr))

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
flog.info("Running with parameters: %s", params, capture = TRUE)

# these options are read from the yml
default_options <- list(
  # The path to the directory containing the cellranger
  # raw matrix and associates barcodes and features
  "cellrangerDir" = "",

  # Path to the emptyDrops output file
  "outdir" = "library.dir"
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

sce <- read10xCounts(opt$cellrangerDir)
flog.info("Dimensions of input object: %s", dim(sce))

sce <- sce[,colSums(counts(sce))>0]
sce <- sce[rowSums(counts(sce))>0,]
flog.info("Dimensions of filtered object: %s", dim(sce))

mean_UMI_per_cell <- data.frame(BARCODE = colData(sce)$Barcode,
                                mean_umi = colMeans(counts(sce)),
                                stringsAsFactors = F)

write.table(mean_UMI_per_cell, file = gzfile(file.path(opt$outdir, "meanUmiPerCell.tsv.gz")),
            sep="\t", quote = FALSE, row.names = FALSE)

flog.info("Finished writing output...")

flog.info("Completed")
