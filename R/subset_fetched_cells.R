# Libraries --------------------------------------------------------------------

stopifnot(require(optparse),
          require(futile.logger),
          require(R.utils),
          require(dplyr),
          require(data.table),
          require(Matrix),
          require(BiocParallel)
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--celltable"), default="fetch.cells.dir/cell.table.tsv.gz",
              help="cell table output of the fetch_cells pipelines."),
  make_option(c("--mincell"), default=100,
              help="Minimum number cells per sample."),
  make_option(c("--samplevar"), default="sample_id",
              help="Sample/batch variable."),
  make_option(c("--exclude"), default="none",
              help="Single column file with barcodes to exclude."),
  make_option(c("--log_filename"), default="subset_fetched_cells.log")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger -----------------------------------------------------------------------
# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -
flog.threshold(INFO)

# now set to append mode -------------------------------------------------------
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

multicoreParam <- MulticoreParam(workers = opt$numcores)

# Read in data -----------------------------------------------------------------

flog.info("Reading cell table ...")
celltable <- fread(file.path(opt$celltable))
flog.info("Number of barcodes in input: %s", format(nrow(celltable),
                                                    big.mark=","))

# Remove unwanted barcodes -----------------------------------------------------

if(opt$exclude != "none" & file.exists(opt$exclude)) {
  
  flog.info("Cell table before excluding barcodes:", 
            format(nrow(celltable), big.mark=","), capture = TRUE)
  
  exclude_barcodes <- fread(opt$exclude, h = F)[[1]]
  
  flog.info("Barcodes to exclude:", 
            format(length(exclude_barcodes), big.mark=","), capture = TRUE)
  
  celltable %>%
    filter(!barcode_id %in% exclude_barcodes) -> celltable
  
  flog.info("Cell table after excluding barcodes:", 
            format(nrow(celltable), big.mark=","), capture = TRUE)
  
}

# Sample subset ----------------------------------------------------------------

if(!opt$samplevar %in% colnames(celltable)) {
  stop("sample_var not present in celltable, please reconsider this parameter")
}

celltable %>%
  group_by(!!rlang::sym(opt$samplevar)) %>%
  summarise('ncell' = n()) -> sample_meta

filter(sample_meta, ncell > opt$mincell)[[opt$samplevar]] -> keep_samples

flog.info("Samples to keep:", keep_samples, capture = TRUE)

# Filter cell table ------------------------------------------------------------

celltable %>%
  filter(!!rlang::sym(opt$samplevar) %in% keep_samples) -> sub_celltable

flog.info("Sub-sampled cell table barcodes:", format(nrow(sub_celltable), 
                                                     big.mark=","), 
          capture = TRUE)

# Write results ----------------------------------------------------------------

outfile <- gsub("\\.gz", "", opt$celltable)

write.table(sub_celltable, outfile, 
            sep = "\t", quote = FALSE, row.names = FALSE)

system(paste("rm", opt$celltable))
gzip(outfile)

flog.info("Sub cell_table has been written.")
