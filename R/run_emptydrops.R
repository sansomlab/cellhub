## Title ----
##
## Run emptydrops to detect empty droplets
##
## Description ----
## Run emptyDrops to assess if droplets are empty or might contain cells
##
## Default parameters ----
## example yml: /gfs/devel/kjansen/dropflow/Rmd/emptydrops.yml

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(optparse))
stopifnot(require(DropletUtils))
stopifnot(require(ggplot2))
stopifnot(require(futile.logger))
stopifnot(require(knitr))
stopifnot(require(optparse))

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

  # FDR filter to use for EmptyDrops
  "FDR" = 0.001,

  # Path to the emptyDrops output file
  "outdir" = "sample.dir"

  # Path to the file containing barcodes for blacklisting
  #"blacklist" = NULL
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


## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

sce <- read10xCounts(file.path(opt$cellrangerDir))
flog.info("Dimensions of input object: %s", dim(sce))

## ######################################################################### ##
## ############## (ii) Filter data including blacklisting ################## ##
## ######################################################################### ##

# if (!is.null(opt$blacklist)) {
#   stopifnot(file.exists(opt$blacklist))
#   flog.info("Read blacklisted cells from: %s", opt$blacklist)
#   cells_blacklist <- data.frame(barcodes_in = scan(opt$blacklist, "character"))
#   cells_blacklist$barcodes_blacklist = paste(cells_blacklist$barcodes_in, "1", sep="-")
#   flog.info("Number of cells in given blacklist: %s", nrow(cells_blacklist))
#   cells_remove <- colData(sce)$Barcode %in% cells_blacklist$barcodes_blacklist
#   sce <- sce[, !cells_remove]
#   flog.info("Number of cells after blacklisting: %s", ncol(sce))
# } else {
#   flog.info("No cells were blacklisted.") }

sce <- sce[,colSums(counts(sce))>0]
sce <- sce[rowSums(counts(sce))>0,]
flog.info("Dimensions of filtered object: %s", dim(sce))

my.counts <- counts(sce)

## ######################################################################### ##
## ################## (iii) Run emptyDrops ################################# ##
## ######################################################################### ##

set.seed(100)
e.out <- emptyDrops(my.counts)

is.cell <- e.out$FDR <= opt$FDR
flog.info("Number of non-empty droplets: %s", sum(is.cell, na.rm=TRUE))

## make input file for rank plots
br.out <- barcodeRanks(counts(sce))
knee_point <- metadata(br.out)$knee
inflection_point <- metadata(br.out)$inflection

br.out <- as.data.frame(br.out)
br.out$barcode <- colData(sce)$Barcode


## ######################################################################### ##
## ###################### (iv) Store output ################################ ##
## ######################################################################### ##

flog.info("Writing output...")

flog.info("Write out ranking table...")
write.table(br.out, file = gzfile(file.path(opt$outdir, "barcode_ranks.tsv.gz")), sep="\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(knee_point = knee_point, inflection_point = inflection_point),
            file = gzfile(file.path(opt$outdir, "barcode_ranks_metadata.tsv.gz")), sep="\t", quote = FALSE, row.names = FALSE)

emptydrops_out = as.data.frame(e.out)
emptydrops_out$barcode = colData(sce)$Barcode
emptydrops_out$emptydrops_cell <- FALSE
emptydrops_out$emptydrops_cell[emptydrops_out$FDR <= opt$FDR] <- TRUE

write.table(emptydrops_out, file = gzfile(file.path(opt$outdir, "emptyDrops.tsv.gz")), sep="\t", quote = FALSE, row.names = FALSE)

flog.info("Finished writing output...")

flog.info("Completed")
