# ------------------------------------------------------------------------------
# DSB normalized ADT UMAP embedding & visualization 
# The input for this task is the fetch_cell pipeline output
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
          require(harmony),
          require(ggplot2),
          require(cowplot),
          require(matrixStats),
          require(ggridges)
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--adt_dsb"), default="fetch.cells.dir/ADT.mtx.dsb.full.dir",
              help="Folder with the filtered ADT DSB nrmalized market matrices."),
  make_option(c("--nfeat"), default="all",
              help="If not all, the number of top most variable proteins."),
  make_option(c("--cellmeta"), default="fetch.cells.dir/cell.table.tsv.gz",
              help="Cell indexed metadata table."),
  make_option(c("--colorvar"), default="library_id",
              help="Variable name to color cells."),
  make_option(c("--numcores"), default=8,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="dsb_umap_embedding.log"),
  make_option(c("--outfile"), default="fetch.cells.dir/ADT.mtx.full.dir/dsb_umap_embedding.pdf")
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

plot_adt_feat <- function(gene, adt, umap_clr, split=FALSE) {
  print(gene)
  clr_exp <- data.frame(adt[gene, ])
  clr_exp[["barcode_id"]] <- rownames(clr_exp)
  colnames(clr_exp)[1] <- gene
  dim(clr_exp)
  head(clr_exp)
  umap_dat <- merge(umap_clr, clr_exp, by = "barcode_id")
  umap_dat <- umap_dat[sample(1:nrow(umap_clr), nrow(umap_clr)*0.2), ]
  print(dim(umap_dat))
  
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
  
  if(split) {
    nc <- ceiling(sqrt(length(unique(umap_dat[[opt$colorvar]]))))
    return(gg_umap + facet_wrap(as.formula(paste0("~", opt$colorvar)), ncol = nc))
  } else {
    return(gg_umap)
  }
}

# Read in data -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Read DSB normalized ADT data -------------------------------------------------
# ------------------------------------------------------------------------------

flog.info("Reading DSB data object...")

if(!dir.exists(opt$adt_dsb)) {
  stop("No folder with market matrices exist.")
}

dsb_adt_mtx <- Read10X(file.path(opt$adt_dsb))

flog.info("Number of barcodes in input: %s", format(ncol(dsb_adt_mtx), 
                                                    big.mark=","))

print(dim(dsb_adt_mtx))

# DSB - UMAP -------------------------------------------------------------------
# umap straight from DSB normalized values -------------------------------------

print(dim(dsb_adt_mtx))
print(head(t(data.matrix(dsb_adt_mtx))))

# Feature selection ------------------------------------------------------------

dsbmtx <- data.matrix(dsb_adt_mtx)

prots <- data.frame('var' = rowVars(dsbmtx),
                    'mean' = rowMeans(dsbmtx),
                    'sum' = rowSums(dsbmtx),
                    'pct' = rowSums(dsbmtx > 0)/ncol(dsbmtx)) %>%
  arrange(desc(var), desc(pct), desc(mean))

prots[["protein"]] <- rownames(prots)

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

pdf(gsub("\\.pdf", paste0("_dsb_protein_variation.pdf"), opt$outfile),
    height = 4, width = 8)

  plot_grid(plotlist = list(gg_feat, gg_feat_rank), nrow = 1)

dev.off()

if(!opt$nfeat == "all") {
  genes <- head(prots, opt$nfeat)[["protein"]]
} else {
  genes <- prots$protein
}

umap <- RunUMAP(t(data.matrix(dsb_adt_mtx[genes, ])),
                dims = 1:length(genes), 
                n.components = 2, 
                min.dist = 0.2)

umap_clr <- data.frame(umap@cell.embeddings)
umap_clr[["barcode_id"]] <- rownames(umap_clr)
head(umap_clr)
dim(umap_clr)

write.table(umap_clr, gsub("\\.pdf", paste0("_DSB_umap_barcode_embedding.tsv"), 
                           opt$outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)

if(!file.exists(opt$cellmeta)) {
  stop("Cell metadata file non existant.")
}

meta <- fread(opt$cellmeta)
print(head(meta))
umap_clr <- merge(umap_clr, meta, by = "barcode_id")
rownames(umap_clr) <- umap_clr$barcode_id

gg_umap <- ggplot(umap_clr, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = !!rlang::sym(opt$colorvar)), 
             size = 0.01, lwd = 0.01) +
  geom_density_2d(lwd = 0.2, color = "black", alpha = .6) +
  theme_void()

png(gsub("\\.pdf", "_umap_dsb_protein.png", opt$outfile),
    height = 2800, width = 2800, res = 300)

  plot(gg_umap)

dev.off()

ggs <- lapply(prots$protein, function(gene) {
  plot_adt_feat(gene, dsb_adt_mtx, umap_clr)
})

nc <- ceiling(sqrt(length(ggs)))

png(gsub("\\.pdf", paste0("_umap_dsb_protein_expression.png"), opt$outfile),
    height = 1000*nc, width = 1000*nc, res = 300)

  plot_grid(plotlist = ggs, nrow = nc)

dev.off()


# UMAP plot spliting by colorvar

lapply(prots$protein, function(gene) {
  print(gene)
  
  png(gsub("\\.pdf", paste0("_umap_dsb_protein_expression_split_", gene, ".png"), opt$outfile),
      height = 1000, width = 1000, res = 300)
    plot(
      plot_adt_feat(gene, dsb_adt_mtx, umap_clr, split=TRUE)
    )
  dev.off()
  
})


# Density plots ----------------------------------------------------------------

sample_ids <- unique(meta[[opt$colorvar]])

gg_dsb <- lapply(sample_ids, function(s) {
  b <- filter(meta, !!rlang::sym(opt$colorvar) == s)[["barcode_id"]]
  sub_dsb <- data.matrix(dsb_adt_mtx[, colnames(dsb_adt_mtx) %in% b])
  dsb <- data.frame(melt(sub_dsb))
  ggplot(dsb, aes(x = value, y = Var1, fill = Var1)) +
    geom_density_ridges_gradient(scale = 3, 
                                 rel_min_height = 0.01, 
                                 gradient_lwd = .01,
                                 alpha = .3) +
    theme_ridges(font_size = 4, grid = TRUE) + 
    theme(axis.title.y = element_blank()) +
    xlim(NA, quantile(dsb$value, 0.98)) +
    xlab("DSB normalized value") +
    theme(legend.position = "none") +
    ggtitle(s)
})

nc <- ceiling(sqrt(length(gg_dsb)))

pdf(gsub(".pdf", paste0("_DSB_density_", opt$colorvar, ".pdf"), opt$outfile), 
    width = 4*nc, height = 0.5*nrow(dsb_adt_mtx)*nc)

  plot_grid(plotlist = gg_dsb, ncol = nc)

dev.off()

# heatmap per colorvar ---------------------------------------------------------

dsb_mat <- data.frame(t(data.matrix(dsb_adt_mtx)))
dsb_mat[["barcode_id"]] <- rownames(dsb_mat)
dsb_mat <- merge(dsb_mat, select(meta, barcode_id, !!rlang::sym(opt$colorvar)),
                 by = "barcode_id")
head(dsb_mat)

dsb_ext <- melt(dsb_mat)
head(dsb_ext)

dsb_ext %>%
  dplyr::group_by(!!rlang::sym(opt$colorvar), variable) %>%
  dplyr::summarise('mean' = mean(value)) -> mean_dsb

head(mean_dsb)  

dsb_cl <- data.frame(dcast(data = mean_dsb, as.formula(paste0(opt$colorvar, 
                                                              "~variable"))))
dsb_cl[1:3, 1:3]

mat <- data.matrix(dsb_cl[, 2:ncol(dsb_cl)])
rownames(mat) <- dsb_cl[[opt$colorvar]]
print(mat[1:3, 1:3])

pdf(gsub(".pdf", paste0("_DSB_heatmap_", opt$colorvar, ".pdf"), opt$outfile),
    width = 2+(0.2*nrow(mat)), height = 0.3*ncol(mat))

pheatmap::pheatmap(t(mat), 
                   scale = "row",
                   fontsize_row = 6, 
                   border_color = NA,
                   na.rm = TRUE,
                   cutree_rows = NA,
                   cutree_cols = NA,
                   cluster_cols = TRUE
)
pheatmap::pheatmap(t(mat[, genes]), 
                   scale = "row",
                   fontsize_row = 6, 
                   border_color = NA,
                   na.rm = TRUE,
                   cutree_rows = NA,
                   cutree_cols = NA,
                   cluster_cols = TRUE
)
pheatmap::pheatmap(t(mat), 
                   scale = "none",
                   fontsize_row = 6, 
                   border_color = NA)

dev.off()