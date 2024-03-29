# Libraries --------------------------------------------------------------------

stopifnot(require(optparse),
          require(futile.logger),
          require(R.utils),
          require(dplyr),
          require(data.table),
          require(Matrix),
          require(BiocParallel),
          require(Seurat),
          require(dsb),
          require(ggplot2),
          require(ggridges)
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--unfiltered_dir"), default=".",
              help="Folder with unfiltered cellranger output. Must include
              barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--filtered_dir"), default=".",
              help="Folder with filtered cellranger output. Must include 
              barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--library_id"), default=NULL,
              help="library or channel id"),
  make_option(c("--rm_feat"), default="None",
              help="Comma separated string with features to remove from
              analysis."),
  make_option(c("--qc_bar"), default="None",
              help="Single column .tsv.gz file with the cell-barcodes with high
              GEX-based QC metrics. (Output fetch_cells pipeline.)"),
  make_option(c("--gex_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--adt_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--cap_thr"), default=0.02,
              help="Top fraction of most & last cell-barcodes to remove from 
              analysis. (Number of genes, Number of UMI.)"),
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

# Logger -----------------------------------------------------------------------
# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -
flog.threshold(INFO)

# now set to append mode -------------------------------------------------------
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

multicoreParam <- MulticoreParam(workers = opt$numcores)

# Read in data -----------------------------------------------------------------

flog.info("Reading unfiltered data object...")
raw_adt_mtx <- Read10X(file.path(opt$unfiltered_dir))
flog.info("Number of barcodes in input: %s", format(dim(raw_adt_mtx), 
                                                    big.mark=","))

if(opt$rm_feat != "None") {
  rmft <- unlist(strsplit(opt$rm_feat, ","))
  if(any(!rmft %in% rownames(raw_adt_mtx))) {
    stop("Features provided to be removed are not present.")
  }
  raw_adt_mtx <- raw_adt_mtx[!rownames(raw_adt_mtx) %in% rmft, ]
}

flog.info("Number of barcodes in input: %s", format(dim(raw_adt_mtx), 
                                                    big.mark=","))

flog.info("Reading filtered data object...")
filt_adt_mtx <- Read10X(file.path(opt$filtered_dir))
flog.info("Number of barcodes in input: %s", format(dim(filt_adt_mtx), 
                                                    big.mark=","))

if(opt$rm_feat != "None") {
  rmft <- unlist(strsplit(opt$rm_feat, ","))
  if(any(!rmft %in% rownames(filt_adt_mtx))) {
    stop("Features provided to be removed are not present.")
  }
  filt_adt_mtx <- filt_adt_mtx[!rownames(filt_adt_mtx) %in% rmft, ]
}

flog.info("Number of barcodes in input: %s", format(dim(filt_adt_mtx), 
                                                    big.mark=","))

# Reading GEX seq-depth metrics ------------------------------------------------
if(!file.exists(opt$gex_depth)) {
	stop("Non gex depth metrics")
}

gex_depth <- fread(opt$gex_depth)
flog.info("Head gex-depth QC table:", head(gex_depth), capture = TRUE) 
flog.info("Groups gex-depth QC table:", table(gex_depth$group), capture = TRUE)

if(opt$qc_bar != "None") {
  qcbar <- fread(opt$qc_bar, h = F)[[1]]
  if(length(which(qcbar %in% gex_depth$barcode_id)) == 0) {
    stop("QC barcodes are not present.")
  }
  gex_depth <- filter(gex_depth, `group` == "background" | barcode_id %in% qcbar)
}
flog.info("Head gex-depth QC table:", head(gex_depth), capture = TRUE) 
flog.info("Groups gex-depth QC table:", table(gex_depth$group), capture = TRUE)

# Reading ADT seq-depth metrics  -----------------------------------------------
if(!file.exists(opt$adt_depth)) {
	stop("Non adt depth metrics")
}

adt_depth <- fread(opt$adt_depth)
flog.info("Head adt-depth QC table:", head(adt_depth), capture = TRUE) 
flog.info("Groups adt-depth QC table:", table(adt_depth$group), capture = TRUE)

if(opt$qc_bar != "None") {
  qcbar <- fread(opt$qc_bar, h = F)[[1]]
  if(length(which(qcbar %in% adt_depth$barcode_id)) == 0) {
    stop("QC barcodes are not present.")
  }
  adt_depth <- filter(adt_depth, `group` == "background" | barcode_id %in% qcbar)
}

flog.info("Head adt-depth QC table:", head(adt_depth), capture = TRUE) 
flog.info("Groups adt-depth QC table:", table(adt_depth$group), capture = TRUE)


# Split unfiltered ADT count matrix into bacground & cell containing matrices --

if(opt$bcmin == "None") {
  
  flog.info("Defining background and cell containing barcodes based on \
            authomated procedure.")
  
  cell_barcodes <- colnames(filt_adt_mtx)

  flog.info("Number of cell barcodes:", length(cell_barcodes), 
            capture = TRUE)
  
  
  # Filter background ----------------------------------------------------------
  
  sub_gex_depth <- gex_depth[, c(1, 6, 7)]
  colnames(sub_gex_depth)[2] <- "log_nfeat_gex"
  head(sub_gex_depth)
  
  subadt_depth <- adt_depth[, c(1, 5)]
  colnames(subadt_depth)[2] <- "log_total_UMI_adt"
  head(subadt_depth)
  
  depth <- merge(sub_gex_depth, 
                 subadt_depth, 
                 by = "barcode_id")
  head(depth)
  flog.info("Dim depth:", dim(depth), capture = TRUE)
  
  depth %>%
    filter(log_total_UMI_adt >= 0.5) -> sub_depth
    #filter(log_nfeat_gex >= 1) 
  
  flog.info("Dims sub depth:", dim(sub_depth), 
            capture = TRUE)
  
  back <- filter(sub_depth, group == "background")
  
  if(nrow(back) < 100) {
    
    flog.info("WARNING: Number of background barcodes with more than
            10 adt UMIs is less than 100.")
    
    depth %>%
      filter(log_total_UMI_adt >= 0.2) -> sub_depth
      #filter(log_nfeat_gex >= 0.2) 
    dim(sub_depth)
    
    back <- filter(sub_depth, group == "background")
    
  }
  
  back$log_total_UMI_adt[is.infinite(back$log_total_UMI_adt)] <- 0
  max_adt <- mean(back$log_total_UMI_adt)+(2*sd(back$log_total_UMI_adt))
  quantile(filter(adt_depth, `group` == "cell")[["log_total_UMI"]], 
           opt$cap_thr, na.rm = TRUE) -> ccmin
  
  if(ccmin == 0) {
    quantile(filter(adt_depth, `group` == "cell")[["log_total_UMI"]], 
             opt$cap_thr + .02, na.rm = TRUE) -> ccmin
  }
  if(max_adt >= ccmin) {
    max_adt <- ccmin-0.2
  }
  
  min_adt <- mean(back$log_total_UMI_adt)-(3*sd(back$log_total_UMI_adt))
  if(min_adt > 0.5) {
    min_adt <- 0.5
  } else if(min_adt < 0) {
    min_adt <- 0.5
  }
  
  back$log_nfeat_gex[is.infinite(back$log_nfeat_gex)] <- 0
  head(back)
  max_gen <- mean(back$log_nfeat_gex)+(2*sd(back$log_nfeat_gex))
  
  if(max_gen > 80) {
    max_gen <- log10(80)
  }
  
  quantile(filter(gex_depth, `group` == "cell")[["nfeat"]], opt$cap_thr, 
           na.rm = TRUE) -> cfmin
  
  if(cfmin < 200) {
    flog.info("WARNING!: Minimum number of genes in cell is less than 200", 
              cfmin, capture = TRUE)
  } 
  
  min_gen <- mean(back$log_nfeat_gex)-(3*sd(back$log_nfeat_gex))
  if(min_gen > 0) {
    min_gen <- 0
  } else if(min_gen < 0) {
    min_gen <- 0
  }
  
  # Filtering out lowly & highly seq cell-barcodes -----------------------------
  quantile(filter(gex_depth, `group` == "cell")[["nfeat"]], 1-opt$cap_thr,
           na.rm = TRUE) -> cfmax
  
  cell <- filter(sub_depth, group == "cell") %>%
    filter(log_total_UMI_adt >= max_adt) %>%
    filter(log_nfeat_gex >= log10(200)) %>%
    filter(log_nfeat_gex <= log10(cfmax))
  
  cell_barcodes <- cell_barcodes[cell_barcodes%in%cell[["barcode_id"]]]
  
  # Filtering out lowly & highly seq background-barcodes -----------------------
  back %>%
    filter(log_total_UMI_adt > min_adt &
             log_total_UMI_adt < max_adt &
             log_nfeat_gex > min_gen &
             log_nfeat_gex < max_gen) -> sub_back
  
  back_barcodes <- sub_back[["barcode_id"]]
  
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

flog.info("Subseting unfiltered adt mtx ...")
back_adt_mtx <- raw_adt_mtx[, back_barcodes]
flog.info("background mat:", dim(back_adt_mtx), capture = T)
cell_adt_mtx <- raw_adt_mtx[, cell_barcodes]
flog.info("cell mat:", dim(cell_adt_mtx), capture = T)

# DSB normalization ------------------------------------------------------------

flog.info("DSBNormalizeProtein ...")

adt_norm <- tryCatch({
  
  DSBNormalizeProtein(
    cell_protein_matrix = cell_adt_mtx, 
    empty_drop_matrix = back_adt_mtx, 
    denoise.counts = FALSE)
  
}, error=function(e) {
  
  flog.info("DSB normalization failed based on the input parameters. 
            Consider re-setting the cell & background definitions")
  flog.info("DSB error:", e, capture = TRUE)
  return(NULL)
  
})

# Write market matrices --------------------------------------------------------

if(!is.null(adt_norm)) {
  
  norm_adt_ext <- melt(adt_norm)
  print(head(norm_adt_ext))
  
  gg_dsb <- ggplot(data.frame(norm_adt_ext),
                      aes(x = value, y = Var1, fill = Var1)) +
    geom_density_ridges_gradient(scale = 3, 
                                 rel_min_height = 0.01, 
                                 gradient_lwd = .01,
                                 alpha = .3) +
    theme_ridges(font_size = 8, grid = TRUE) + 
    theme(axis.title.y = element_blank()) +
    xlab("DSB normalized expression") +
    theme(legend.position = "none") +
    ggtitle("DSB normalization")
  
  pdf(gsub("\\.mtx", "_dsb_norm_density.pdf", opt$outfile),
      width = 6, height = 0.15*nrow(adt_norm))
  
  plot(
    gg_dsb
  )
  
  dev.off()
  
  flog.info("Writing output table")
  
  barcodes <- data.frame('f' = colnames(adt_norm))
  bfile <- paste0(dirname(opt$outfile), '/barcodes.tsv')
  write.table(barcodes, bfile, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  gzip(bfile)
  
  features <- data.frame('Id' = rownames(adt_norm))
  features[['Name']] <- rownames(adt_norm)
  features[["Class"]] <- "Antibody Capture"
  ffile <- paste0(dirname(opt$outfile), '/features.tsv')
  write.table(features, ffile, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  gzip(ffile)
  
  sparse_adt_norm_tab <- Matrix(adt_norm, sparse = TRUE)
  writeMM(obj = sparse_adt_norm_tab, file = opt$outfile)
  gzip(opt$outfile)
  
  flog.info("ADT normalized matrix written.")
  
} else {
  
  flog.info("ADT normalization failed for this sample.")
  
}

