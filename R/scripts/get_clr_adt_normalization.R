# ------------------------------------------------------------------------------
# Get CLR normalization ADT 
# The input for this task is the fetch_cell pipeline output
# Save martket matrices with CLR normalized values
# ------------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse),
          require(futile.logger),
          require(R.utils),
          require(dplyr),
          require(data.table),
          require(Matrix),
          require(BiocParallel),
          require(Seurat),
          require(ggplot2),
          require(cowplot),
          require(matrixStats),
          require(ggridges)
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--adt"), default="ADT.mtx.full.dir",
              help="Folder with the filtered ADT DSB nrmalized market matrices."),
  make_option(c("--nfeat"), default="all",
              help="If not all, the number of top most variable proteins."),
  make_option(c("--rm_feat"), default="None",
              help="Comma separated string with features to remove from
              analysis."),
  make_option(c("--qc_bar"), default="None",
              help="Single column .tsv.gz file with the cell-barcodes with high
              GEX-based QC metrics. (Output fetch_cells pipeline.)"),
  make_option(c("--numcores"), default=8,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="clr_adt_norm.log"),
  make_option(c("--outfile"), default="matrix.mtx")
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

# Aux function -----------------------------------------------------------------

plot_adt_feat <- function(gene, adt, umap_clr) {
  print(gene)
  clr_exp <- data.frame(adt@assays$RNA@data[gene, ])
  clr_exp[["barcode_id"]] <- rownames(clr_exp)
  colnames(clr_exp)[1] <- gene
  dim(clr_exp)
  head(clr_exp)
  umap_dat <- merge(umap_clr, clr_exp, by = "barcode_id")
  head(umap_dat)
  
  gg_umap <- ggplot(umap_dat, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = !!rlang::sym(gene)), 
               size = 0.01, lwd = 0.01) +
    theme_void() +
    theme(legend.position = "bottom") +
    scale_colour_gradientn(colors = c("grey95",
                                      "#49006a"), 
                           limits = c(0, quantile(umap_dat[[gene]], .98)),
                           na.value = "grey90") +
    ggtitle(gene)
  
  return(gg_umap)
}


# Read in data -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Read UMI count ADT data ------------------------------------------------------
# ------------------------------------------------------------------------------

flog.info("Reading ADT counts data object...")

if(!dir.exists(opt$adt)) {
  stop("No folder with market matrices exist.")
}

adt_mtx <- Read10X(file.path(opt$adt))

flog.info("Number of barcodes in input: %s", format(ncol(adt_mtx), 
                                                    big.mark=","))
flog.info("Number of features in input: %s", format(nrow(adt_mtx), 
                                                    big.mark=","))

if(opt$rm_feat != "None") {
  rmft <- unlist(strsplit(opt$rm_feat, ","))
  if(any(!rmft %in% rownames(adt_mtx))) {
    stop("Features provided to be removed are not present.")
  }
  adt_mtx <- adt_mtx[!rownames(adt_mtx) %in% rmft, ]
  flog.info("Number of features in input: %s", format(dim(adt_mtx), 
                                                      big.mark=","))
}

if(opt$qc_bar != "None") {
  qcbar <- fread(opt$qc_bar, h = F)[[1]]
  if(length(which(qcbar %in% colnames(adt_mtx))) == 0) {
    flog.info("No QC barcodes present in this sample.")
    quit("no")
  }
  adt_mtx <- adt_mtx[, colnames(adt_mtx) %in% qcbar]
  flog.info("Number of QC barcodes in input: %s", format(ncol(adt_mtx), 
                                                      big.mark=","))
  
}

# Feature selection ------------------------------------------------------------

dsbmtx <- data.matrix(adt_mtx)

prots <- data.frame('var' = rowVars(dsbmtx),
                    'mean' = rowMeans(dsbmtx),
                    'sum' = rowSums(dsbmtx),
                    'pct' = rowSums(dsbmtx > 0)/ncol(dsbmtx)) %>%
  arrange(desc(var), desc(pct), desc(mean))

prots[["protein"]] <- rownames(adt_mtx)

flog.info("Protein/feature stats: %s", format(head(prots), big.mark=","))

gg_feat <- ggplot(prots, aes(x = `mean`, 
                             y = `var`, 
                             color = pct,
                             label = rownames(prots))) +
  geom_point(alpha = 0.5) +
  scale_color_gradient2(low = "grey80", high = "darkred") +
  geom_text(check_overlap = TRUE, size = 1, nudge_x = 1, hjust = 0) +
  theme_classic()

gg_feat_rank <- ggplot(prots, aes(x = `protein`,
                                  y = `var`, 
                                  color = pct)) +
  geom_point(alpha = 0.5) +
  scale_x_discrete(limits = prots$protein) +
  scale_color_gradient2(low = "grey80", high = "darkred") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, vjust=0.5, 
                                   hjust=1, size = 4))

pdf(gsub("\\.mtx", paste0("_UMI_protein_variation.pdf"), opt$outfile),
    height = 4, width = 8)

  plot_grid(plotlist = list(gg_feat, gg_feat_rank), nrow = 1)

dev.off()

if(opt$nfeat != "all") {
  genes <- head(prots, opt$nfeat)[["protein"]]
} else {
  genes <- prots$protein
}

# -- CLR normalization ---------------------------------------------------------

flog.info("Seurat object creation...")

adt <- CreateSeuratObject(counts = adt_mtx[genes, ])
adt <- NormalizeData(adt,
                     normalization.method = 'CLR', 
                     margin = 2)

norm_adt <- data.matrix(adt@assays$RNA@data)

flog.info(norm_adt[1:3, 1:4])

# Write market matrices --------------------------------------------------------

if(!is.null(norm_adt)) {
  
  flog.info("Writing output table")
  
  barcodes <- data.frame('f' = colnames(norm_adt))
  bfile <- paste0(dirname(opt$outfile), '/barcodes.tsv')
  write.table(barcodes, bfile, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  gzip(bfile)
  
  features <- data.frame('Id' = rownames(norm_adt))
  features[['Name']] <- rownames(norm_adt)
  features[["Class"]] <- "Antibody Capture"
  ffile <- paste0(dirname(opt$outfile), '/features.tsv')
  write.table(features, ffile, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  gzip(ffile)
  
  sparse_adt_norm_tab <- Matrix(norm_adt, sparse = TRUE)
  writeMM(obj = sparse_adt_norm_tab, file = opt$outfile)
  gzip(opt$outfile)
  
  flog.info("ADT normalized matrix written.")
  
} else {
  
  flog.info("ADT normalization failed for this sample.")
  
}

norm_adt_ext <- melt(norm_adt)
print(head(norm_adt_ext))

gg_dsb <- ggplot(data.frame(norm_adt_ext),
                 aes(x = value, y = Var1, fill = Var1)) +
  geom_density_ridges_gradient(scale = 3, 
                               rel_min_height = 0.01, 
                               gradient_lwd = .01,
                               alpha = .3) +
  theme_ridges(font_size = 8, grid = TRUE) + 
  theme(axis.title.y = element_blank()) +
  xlab("CLR normalized expression") +
  theme(legend.position = "none") +
  ggtitle("CLR normalization")

pdf(gsub("\\.mtx", "_CLR_norm_density.pdf", opt$outfile),
    width = 6, height = 0.15*nrow(norm_adt))

plot(
  gg_dsb
)

dev.off()

