# Librarie ------
stopifnot(require(optparse),
          require(futile.logger),
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
              help="ADT mininum number of detected features of the background 
              interval"),
  make_option(c("--bfmax"), default=NULL,
              help="ADT maxinum number of detected features of the background 
              interval"),
  make_option(c("--ccmin"), default=NULL,
              help="ADT UMI count mininum limit of the cell-containing 
              interval"),
  make_option(c("--ccmax"), default=NULL,
              help="ADT UMI count maxinum limit of the cell-containing 
              interval"),
  make_option(c("--cfmin"), default=NULL,
              help="ADT mininum number of detected features of the 
              cell-containing interval"),
  make_option(c("--cfmax"), default=NULL,
              help="ADT maxinum number of detected features of the 
              cell-containing interval"),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="dsb_norm.log"),
  make_option(c("--outfile"), default="./matrix.mxt.gz")
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

if(is.null(opt$bcmin)) {
  
  gex_cell_barcodes <- filter(gex_depth, group == "cell")[["barcode_id"]]
  adt_cell_barcodes <- filter(adt_depth, group == "cell")[["barcode_id"]]
  cell_barcodes <- intersect(gex_cell_barcodes, adt_cell_barcodes)
  
  gex_back_barcodes <- filter(gex_depth, group == "background")[["barcode_id"]]
  back_barcodes <- gex_back_barcodes[!gex_back_barcodes %in% cell_barcodes]
  
} else {
  
  gex_cell_barcodes <- filter(gex_depth, group == "cell"))[["barcode_id"]]
  
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
bfile <- paste0(dirname(outfile), '/barcodes.tsv')
write.table(barcodes, bfile, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
gzip(bfile)

features <- data.frame('f' = rownames(adt_norm))
ffile <- paste0(dirname(outfile), '/features.tsv')
write.table(features, ffile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
gzip(ffile)

sparse_adt_norm_tab <- Matrix(adt_norm, sparse = TRUE)
writeMM(obj = sparse_adt_norm_tab, file = outfile)
gzip(outfile)

flog.info("ADT normalized matrix written.")