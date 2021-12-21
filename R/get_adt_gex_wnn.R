# ------------------------------------------------------------------------------
# Get ADT + GEX UMAP embeddings 
# The input for this task is the fetch_cell pipeline DSB ADT output (mtx) and 
# the output of the gex integration pipeline (*.h5ad)
# ------------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------

stopifnot(require(Seurat),
          require(optparse),
          require(futile.logger),
          require(Matrix),
          require(data.table),
          require(dplyr),
          require(ggplot2),
          require(reticulate),
          require(anndata),
          require(cowplot),
          require(matrixStats)
)

# Global options ---------------------------------------------------------------

options(stringsAsFactors = F)

# Script arguments -------------------------------------------------------------

option_list <- list(
  make_option(c("--gex_adata"), default="normalized_integrated_anndata.h5ad",
              help="h5ad GEX data object."),
  make_option(c("--low_dim"), default="harmony",
              help="Low dimensional representation to integrate."),
  make_option(c("--ndim"), default=30,
              help="Number of GEX low-dim components to consider."),
  make_option(c("--adt_dsb"), default="ADT.mtx.full.dir",
              help="Folder with the filtered ADT DSB nrmalized market matrices."),
  make_option(c("--nfeat"), default="all",
              help="If not all, the number of top most variable proteins."),
  make_option(c("--cellmeta"), default="cell.table.tsv.gz",
              help="Cell indexed metadata table."),
  make_option(c("--colorvar"), default="library_id",
              help="Variable name to color cells."),
  make_option(c("--virtualenv"), default="cgat_venv",
              help="Python env with adata module."),
  make_option(c("--python"), default="cgat_venv/bin/python",
              help="Python execution path."),
  make_option(c("--numcores"), default=8,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="adt_gex_co_embedding.log"),
  make_option(c("--outfile"), default="wnn_adt_gex.png")
)

opt <- parse_args(OptionParser(option_list=option_list))

# python config ----------------------------------------------------------------

use_virtualenv(opt$virtualenv)
use_python(opt$python)

# Logger -----------------------------------------------------------------------
# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -

flog.threshold(INFO)

# now set to append mode -------------------------------------------------------

flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

# GEX lower dimensional representation -----------------------------------------

if(!file.exists(opt$gex_adata)) {
  stop("No adata file available")
}
adata <- read_h5ad(opt$gex_adata)

lowdim <- adata$obsm[[paste0('X_', opt$low_dim)]]
flog.info("GEX low-dimensional dim:", dim(lowdim), capture = TRUE)

colnames(lowdim) <- paste0("GEXPC_", 1:ncol(lowdim))
rownames(lowdim) <- adata$obs_names
flog.info("GEX low-dimensional:", lowdim[1:3, 1:3], capture = TRUE)

gex_lowdim <- CreateDimReducObject(embeddings = lowdim, assay = "RNA")

# GEX mock UMI matrix to create Seurat object-----------------------------------

umi <- data.matrix(adata$T)
print(dim(umi))
print(umi[1:3, 1:3])

mca <- CreateSeuratObject(counts = umi[1:2, ], assay = "RNA")
print(mca)

mca[["GEXPC"]] <- gex_lowdim

# ADT DSB normalized matrix ----------------------------------------------------

if(!dir.exists(opt$adt_dsb)) {
  stop("No DSB normalized folder.")
}
dsb <- Read10X(opt$adt_dsb)

# Protein variance filter
dsbmtx <- data.matrix(dsb)
prots <- data.frame('var' = rowVars(dsbmtx),
                    'mean' = rowMeans(dsbmtx),
                    'sum' = rowSums(dsbmtx),
                    'pct' = rowSums(dsbmtx > 0)/ncol(dsbmtx)) %>%
  arrange(desc(var), desc(pct), desc(mean))
prots[["protein"]] <- rownames(prots)

dsb <- t(data.matrix(dsb))
print(dim(dsb))
print(dsb[1:3, 1:3])

if(!opt$nfeat == "all") {
  genes <- head(prots, opt$nfeat)[["protein"]]
} else {
  genes <- prots$protein
}

dsb <- dsb[, genes]

feat_pcs <- data.frame('prot' = colnames(dsb),
                       'name' = paste0("ADTPC_", 1:ncol(dsb)))

print(head(feat_pcs))
write.table(feat_pcs, gsub("\\.png", "_prot2lowdim.tsv", opt$outfile),
            sep = "\t", quote = FALSE, row.names = FALSE)

colnames(dsb) <- feat_pcs$name

adt_lowdim <- CreateDimReducObject(embeddings = dsb,
                                   assay = "ADT")

adt_assay <- CreateAssayObject(data = t(dsb[colnames(umi), ]))

mca[["ADT"]] <- adt_assay

mca[["ADTPC"]] <- adt_lowdim

# run WNN ----------------------------------------------------------------------

mca = FindMultiModalNeighbors(
  object = mca,
  reduction.list = list("GEXPC", "ADTPC"),
  weighted.nn.name = "dsb_wnn",
  knn.graph.name = "dsb_knn",
  modality.weight.name = "dsb_weight",
  snn.graph.name = "dsb_snn",
  dims.list = list(1:opt$ndim, 1:ncol(dsb)),
  verbose = TRUE
)

print(mca)

# Get UMAP representation ------------------------------------------------------

mca = RunUMAP(mca, 
              nn.name = "dsb_wnn", 
              reduction.name = "dsb_wnn_umap",
              reduction.key = "dsbwnnUMAP_", 
              seed.use = 666, 
              verbose = TRUE)

umap <- data.frame(mca@reductions$dsb_wnn_umap@cell.embeddings)
umap$barcode_id <- rownames(umap)
write.table(umap, gsub("\\.png", "_wnn_adt_gex_umap.tsv", opt$outfile),
            sep = "\t", quote = FALSE, row.names = FALSE)

if(!file.exists(opt$cellmeta)) {
  stop("Cell metadata file non existant.")
}

meta <- fread(opt$cellmeta)
print(head(meta))
umap <- merge(umap, meta, by = "barcode_id")
rownames(umap) <- umap$barcode_id

gg_umap <- ggplot(umap, aes(x = dsbwnnUMAP_1, y = dsbwnnUMAP_2)) +
  geom_point(aes(color = !!rlang::sym(opt$colorvar)), 
             size = 0.01, lwd = 0.01) +
  geom_density_2d(lwd = 0.2, color = "black", alpha = .6) +
  theme_void()

png(gsub("\\.png", paste0("_umap_dsb_", opt$colorvar, ".png"), opt$outfile),
    height = 2800, width = 2800, res = 300)

	plot(gg_umap)

dev.off()

saveRDS(mca, gsub("\\.png", "_wnn_adt_gex_umap.rds", opt$outfile))

# Get clusters -----------------------------------------------------------------

get.clusters <- function(mca, res) {
  cls <- Reduce(cbind,
                lapply(res, function(r){
                  print(r)
                  mca <- FindClusters(mca,
                                      graph.name = "dsb_knn", 
                                      algorithm = 3, 
                                      resolution = r,
                                      random.seed = 666,  
                                      verbose = TRUE)
                  
                  cl <- mca@meta.data[, paste0("dsb_knn_res.", r)]
                  cl <- data.frame(cl)
                  colnames(cl) <- paste0("dsb_wnn_res.", r)
                  cl
                })
  )
  cls[["barcode_id"]] <- rownames(mca@meta.data)
  cls
}

cluster_df <- get.clusters(mca, c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2))
print(head(cluster_df))
write.table(cluster_df, gsub("\\.png", "_wnn_adt_gex_clusters.tsv", opt$outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)


umap <- merge(umap, cluster_df, by = "barcode_id")

cls <- colnames(cluster_df)[-1]

ggs <- lapply(cls, function(cl) {
  umap[[cl]] <- factor(umap[[cl]])
  gg_umap <- ggplot(umap, aes(x = dsbwnnUMAP_1, y = dsbwnnUMAP_2)) +
    geom_point(aes(color = !!rlang::sym(cl)), 
               size = 0.01, lwd = 0.01) +
    theme_void() +
    theme(legend.position = "none") +
    scale_color_brewer(palette="Dark2") +
    ggtitle(cl)
  
  return(gg_umap)
})

nc <- ceiling(sqrt(length(ggs)))

png(gsub("\\.png", paste0("_umap_adt_gex_wnn_clustering.png"), opt$outfile),
    height = 1000*nc, width = 1000*nc, res = 300)

  plot_grid(plotlist = ggs, nrow = nc)

dev.off()