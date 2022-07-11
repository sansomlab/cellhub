## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(openxlsx),
  require(optparse),
  require(ComplexHeatmap),
  require(cellhub)
)

# Options ----

option_list <- list(
    make_option(c("--loom"), default="data.loom",
                help="path to the loom file"),
    make_option(c("--matrix_loc"), default="matrix",
                help="location of the matrix in the loom file"),
    make_option(c("--scale"), default=FALSE,
                help="should the expression data be scaled for the heatmap"),
    make_option(c("--barcode_id_loc="), default="col_attrs/barcode_id",
                help="location of the barcode ids in the loom file"),
    make_option(c("--gene_id_loc="), default="row_attrs/gene_ids",
                help="location of the gene ids in the loom file"),
    make_option(c("--marker_files"), default="none",
                help="List of marker files to aggregate"),
    make_option(c("--metadata"), default="none",
                help="A tsv file containing the metadata"),
    make_option(c("--annotation"), default="none",
                help="A tsv file containing the annotation"),
    make_option(c("--clusterids"), default="none",
                help="A tsv file containing the cluster identities"),
    make_option(c("--subgroup"), default=NULL,
                help="Optional. Name of a column in metadata to add annotation for in the heatmap"),
    make_option(c("--pdf"), default = FALSE,
                help="Create a pdf version of the top marker heatmap"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")
if(!is.null(opt$subgroup)) { opt$subgroup <- strsplit(opt$subgroup,",")[[1]]}
cat("Running with options:\n")
print(opt)

readMarkers <- function(marker_file_name, 
                        min_pct=0.1, 
                        min_fc=1.5)
{
  fn <- basename(marker_file_name)
  fn_parts <- strsplit(fn, "\\.")[[1]]
  
  cluster = fn_parts[1]
  level = fn_parts[2]
  
  x <- read.table(marker_file_name,sep="\t", header=TRUE)
  # colnames: gene_id pval    logfoldchange   pct     pct_other      mean    mean_other  FC      min_FC
  
  x$logfoldchange <- NULL
  
  x$cluster <- as.numeric(cluster)
  x$level <- level
  
  # filter the markers by setting pvals to NA
  x$pval[abs(log2(x$FC)) < log2(min_fc)] <- NA
  
  x$pval[x$pct < min_pct & x$pct_other < min_pct] <- NA
  
  # compute within cluster (& level) p.adjust
  x$p.adj <- p.adjust(x$pval, method="BH")
  
  x
}

mfs <- strsplit(opt$marker_files, ",")[[1]] 

cols <- c("gene_id","pval","p.adj",#"logfoldchange",
          "pct","pct_other","mean_exprs", "mean_exprs_other","FC","min_FC",
          "cluster", "level")

results <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(results) <- cols

for(mf in mfs)
{
  message("Reading markers from: ", mf)
  x <- readMarkers(mf,min_pct=0.1,min_fc=1.5)
  x <- x[,cols]
  results <- rbind(results,x)
}

# write.table(results, gzfile("raw.tsv.gz"), col.names=T, sep="\t", quote=F)

message("Summarising markers across levels")
print(head(results))
# collapse levels 
# - this is where conserved markers are identified
#   (maximum p-value per level is retained)
#   (max(NA, float) = NA)
markers <- results %>% group_by(cluster, gene_id) %>%
   summarize(across(c("pval","p.adj"),max),
             across(c("pct","pct_other", "mean_exprs", "mean_exprs_other", "FC","min_FC"), mean),
             nlevels=n()) %>%
   arrange(cluster,pval,desc(FC))

message("Marker aggregation complete")
markers <- data.frame(markers)
print(dim(markers))

# add the gene names
anno <- parseBiomartAnnotation(opt$annotation)
markers$gene_name <- getGeneNames(anno, markers$gene_id)

#markers <- markers[,c("cluster","gene","gene_id",
#                      "p_val","p.adj",
#                      "avg_logFC","pct.1","pct.2",
#                      "cluster_mean","other_mean")]

markers <- markers[order(markers$cluster, markers$p.adj),]

# markers$index <- c(rownames(markers))


markers$log2FC <- log2(markers$FC)
markers$min_log2FC <- log2(markers$min_FC)

## write out the full (unannotated) table of significantly
## differentially expressed marker genes as a.tsv file.

message("Saving marker.summary.table.tsv")

marker_file <- paste(outPrefix,"table","tsv","gz",
                     sep=".")

write.table(markers,
            gzfile(marker_file),
            quote=F,sep="\t",row.names=F)

################# <-------------------------> ##################

print("Filtering by adjusted p-value")
filtered_markers <- markers[markers$p.adj < 0.1 & !is.na(markers$p.adj),]
message("Filtered markers are:")
print(dim(filtered_markers))

## write out a shorter table of significantly differentially
## expressed marker genes in excel format.
## TODO allow specification of adjusted p-value threshold.

message("Saving marker.summary.table.xlsx")
wb <- createWorkbook()

addWorksheet(wb,"filtered_markers")
setColWidths(wb,"filtered_markers",cols=1:ncol(filtered_markers),widths=10)
hs <- createStyle(textDecoration = "BOLD")
writeData(wb, "filtered_markers", tidyNumbers(filtered_markers),
          withFilter = T, headerStyle=hs)
saveWorkbook(wb, file=paste(outPrefix,"table","xlsx",
                            sep="."),
             overwrite=T)

# message("Making a heatmap of the top marker genes from each cluster")
    
# mch <- markerComplexHeatmap(loom_path=opt$loom,
#                             matrix_loc=opt$matrix_loc,
#                             barcode_id_loc=opt$barcode_id_loc, #"col_attrs/barcode_id",
#                             gene_id_loc=opt$gene_id_loc, # "row_attrs/gene_ids",
#                             scale=opt$scale,
#                             cluster_ids=opt$clusterids,
#                             metadata_file=opt$metadata,
#                             marker_table=filtered_markers,
#                             priority="log2FC",
#                             n_markers=20,
#                             cells_use=NULL,
#                             row_names_gp=10,
#                             sub_group=opt$subgroup)

# drawHeatmap <- function()
# {
#     draw(mch)
# }

# # } else {

# #     drawHeatmap <- function()
# #     {
# #     plot.new()
# #     text(0.5,0.5,"scale.data slot not present")
# #     }
# # }



# save_plots(paste(outPrefix,"heatmap", sep="."),
#            plot_fn=drawHeatmap,
#            width = 7,
#            to_pdf = opt$pdf,
#            height = 9)



## summarise the number of marker genes identified for each cluster

# shouldn't be reading this twice!.
cids <- read.table(gzfile(opt$clusterids), header=T, sep="\t")

summary <- c()
for(id in unique(filtered_markers$cluster))
{
    ncells = nrow(cids[cids$cluster_id==id,])
    npos = length(filtered_markers$p.adj[filtered_markers$cluster==id & filtered_markers$avg_logFC > 0] )
    nneg = length(filtered_markers$p.adj[filtered_markers$cluster==id & filtered_markers$avg_logFC < 0] )
    ntotal = npos + nneg
    summary <- c(summary,c(id, ncells, npos, nneg, ntotal))
}

sumdf <- data.frame(matrix(summary,ncol=5,byrow=T))
colnames(sumdf) <-c("cluster","ncells","n_pos_markers","n_neg_markers","n_markers")

print("Writing out some summary statistics")
write.table(sumdf,
            paste(outPrefix,"stats","tsv",
                  sep="."),
            quote=F,sep="\t",row.names=F)


message("completed")
