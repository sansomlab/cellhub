# Libraries ------
stopifnot(require(optparse),
          require(futile.logger),
          require(R.utils),
          require(dplyr),
          require(data.table),
          require(Matrix),
          require(BiocParallel),
          require(Seurat),
          require(dsb)
)

# Global options ------

options(stringsAsFactors = F)

# Script arguments ------

option_list <- list(
  make_option(c("--cellranger_dir"), default=".",
              help="Folder with filtered cellranger output. Must include 
              barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--library_id"), default=NULL,
              help="library or channel id"),
  make_option(c("--gex_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--adt_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--bcmin"), default=NULL,
              help="ADT UMI count mininum limit of the background interval"),
  make_option(c("--bcmax"), default=NULL,
              help="ADT UMI count maximum limit of the background interval"),
  make_option(c("--bfmin"), default=NULL,
              help="GEX mininum number of detected features of the background 
              interval"),
  make_option(c("--bfmax"), default=NULL,
              help="GEX maxinum number of detected features of the background 
              interval"),
  make_option(c("--ccmin"), default=NULL,
              help="ADT UMI count mininum limit of the cell-containing 
              interval"),
  make_option(c("--ccmax"), default=NULL,
              help="ADT UMI count maxinum limit of the cell-containing 
              interval"),
  make_option(c("--cfmin"), default=NULL,
              help="GEX mininum number of detected features of the 
              cell-containing interval"),
  make_option(c("--cfmax"), default=NULL,
              help="GEX maxinum number of detected features of the 
              cell-containing interval"),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="dsb_norm.log"),
  make_option(c("--outfile"), default="./matrix.mxt")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger ------

# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -
###################
flog.threshold(INFO)

# now set to append mode -------------------------------------------------------
###################
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)


multicoreParam <- MulticoreParam(workers = opt$numcores)

# Read in data -----------------------------------------------------------------
###################

flog.info("Reading data object...")
raw_adt_mtx <- Read10X(file.path(opt$cellranger_dir))
flog.info("Number of barcodes in input: %s", format(ncol(raw_adt_mtx), 
                                                    big.mark=","))

gex_depth <- fread(opt$gex_depth)
adt_depth <- fread(opt$adt_depth)

# Split unfiltered ADT count matrix into bacground & cell containing matrices --
###################

if(opt$bcmin == "None") {
  
  flog.info("Defining background and cell containing barcodes based on \
            authomated procedure.")
  
  gex_cell_barcodes <- filter(gex_depth, group == "cell")[["barcode_id"]]
  adt_cell_barcodes <- filter(adt_depth, group == "cell")[["barcode_id"]]
  cell_barcodes <- intersect(gex_cell_barcodes, adt_cell_barcodes)

  flog.info("Number of cell barcodes:", length(cell_barcodes), 
            capture = TRUE)
  
  gex_back_barcodes <- filter(gex_depth, group == "background")[["barcode_id"]]
  back_barcodes <- gex_back_barcodes[!gex_back_barcodes %in% cell_barcodes]
  
  flog.info("Number of background barcodes:", length(back_barcodes), 
            capture = TRUE)
  
} else {
  
  flog.info("Defining background and cell containing barcodes based on \
            user's thresholds")
  
  gex_cell_barcodes <- filter(gex_depth,
                              ngenes > opt$cfmin &
                              ngenes < opt$cfmax)[["barcode_id"]]
  
  adt_cell_barcodes <- filter(adt_depth,
                              log_total_UMI > opt$ccmin &
                              log_total_UMI < opt$ccmax)[["barcode_id"]]
  
  cell_barcodes <- intersect(gex_cell_barcodes, adt_cell_barcodes)
  
  flog.info("Number of cell barcodes:", length(cell_barcodes), 
            capture = TRUE)
  
  adt_back_barcodes <- filter(adt_depth, 
                              log_total_UMI > opt$bcmin &
                              log_total_UMI < opt$bcmax)[["barcode_id"]]
  
  back_barcodes <- adt_back_barcodes[!adt_back_barcodes %in% cell_barcodes]
  
  flog.info("Number of background barcodes:", length(back_barcodes), 
            capture = TRUE)
  
}

back_adt_mtx <- raw_adt_mtx[, back_barcodes]
cell_adt_mtx <- raw_adt_mtx[, cell_barcodes]

# DSB normalization ------------------------------------------------------------
###################

flog.info("DSBNormalizeProtein ...")

adt_norm = DSBNormalizeProtein(
  cell_protein_matrix = back_adt_mtx, 
  empty_drop_matrix = cell_adt_mtx, 
  denoise.counts = TRUE 
)

# Write market matrices --------------------------------------------------------
###################

flog.info("Writing output table")

barcodes <- data.frame('f' = colnames(adt_norm))
bfile <- paste0(dirname(opt$outfile), '/barcodes.tsv')
write.table(barcodes, bfile, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
gzip(bfile)

features <- data.frame('f' = rownames(adt_norm))
ffile <- paste0(dirname(opt$outfile), '/features.tsv')
write.table(features, ffile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
gzip(ffile)

sparse_adt_norm_tab <- Matrix(adt_norm, sparse = TRUE)
writeMM(obj = sparse_adt_norm_tab, file = opt$outfile)
gzip(opt$outfile)

flog.info("ADT normalized matrix written.")