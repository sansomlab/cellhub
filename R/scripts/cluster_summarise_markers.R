## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(dplyr),
  require(matrixStats),
  require(openxlsx),
  require(optparse),
  require(cellhub)
)

# Options ----

option_list <- list(
    make_option(c("--marker_files"), default="none",
                help="List of marker files to aggregate"),
    make_option(c("--minpct"), default=0.1, type="double",
                help="minimum fraction of cells expressing gene"),
    make_option(c("--minfc"), default=1.5, type="double",
                help="minimum foldchange"),
    make_option(c("--clusterids"), default="none",
                help="A tsv file containing the clusterids"),
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

cols <- c("gene_name","gene_id", "score", "pval","p.adj", #"logfoldchange",
          "pct","pct_other","mean_exprs", "mean_exprs_other","FC","min_FC",
          "cluster", "level")

results <- data.frame(matrix(ncol=length(cols),nrow=0))
colnames(results) <- cols

for(mf in mfs)
{
  message("Reading markers from: ", mf)
  x <- readMarkers(mf,
                   min_pct=opt$minpct,
                   min_fc=opt$minfc)
  
  if(!"gene_id" %in% colnames(x))
  {
    # add fake gene_ids if not present to avoid multiple if clauses
    # downstream. For pathway analysis real gene_ids must be supplied
    # in adata.var.gene_id
    x$gene_id <- x$gene_name
   }
  
  x <- x[,cols]
  results <- rbind(results,x)
}

message("Summarising markers across levels")

print(head(results))

#only keep markers where FC is in same direction
markers.dir <- results %>% group_by(cluster, gene_name, gene_id) %>%
                          mutate(same.dir = ifelse(sum(abs(log2(FC)))==abs(sum((log2(FC)))), TRUE, NA )) %>% ungroup() %>%
                          filter(same.dir)

# collapse levels 
# - this is where conserved markers are identified
#   (maximum p-value per level is retained)
#   (max(NA, float) = NA)

markers <- markers.dir %>% group_by(cluster, gene_name, gene_id) %>%
   summarize(across(c("score"),min),
             across(c("pval","p.adj"),max),             
             across(c("pct","pct_other", "mean_exprs", "mean_exprs_other", "FC","min_FC"), mean),
             nlevels=n()
             ) %>%
   arrange(cluster,pval,desc(FC))
   
# When dealing with conserved markers it is possible that some levels
# had no cells so the DE test was not run.

# detect the number of levels expected as the maximum observed.
# TODO: make this an explicit check using the metadata.
n_levels <- max(markers$nlevels)

# nuke the adjusted p-values for genes not tested across
# all the levels of the conserved factor.
markers$p.adj[markers$nlevels < n_levels] <- NA

message("Marker aggregation complete")
markers <- data.frame(markers)
print(dim(markers))

# add the gene names
markers <- markers[order(markers$cluster, markers$p.adj),]

markers$log2FC <- log2(markers$FC)
markers$min_log2FC <- log2(markers$min_FC)

## write out the full table of significantly
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

message("summarising the number of marker genes identified for each cluster")

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

message("preparing the gene universe")

# The gene universe is defined per-cluster as the set of genes
# expressed > opt$min_pct in the "cluster" and/or
# "other" population.
#
# TODO: size of universe can be very small, consider increasing.
#
# When dealing with conserved markers, markers need to pass
# the min_pct in all levels in order to be included in the results.
# Hence we also use the minimum min_pct for inclusion in the universe
# list.
#
universe <- results %>% group_by(cluster, gene_id) %>%
   summarize(across(c("pct","pct_other"), min),
             nlevels=n()) 

# ignore genes not found in all the levels
n_levels <- max(universe$nlevels)
universe <- universe[universe$nlevels == n_levels, ]

# take the max of the cluster or the other cells
universe$max_pct <- rowMaxs(as.matrix(universe[, c("pct", "pct_other")]))

universe <- universe[universe$max_pct >= opt$minpct,]

for(cluster in unique(universe$cluster))
{
  cluster_universe <- universe[universe$cluster==cluster, "gene_id"] 
  
  message("Saving gene universe for cluster ",cluster)
  write.table(cluster_universe,
              gzfile(file.path(opt$outdir,
                               paste(cluster, "universe.tsv.gz", sep="."))),
              quote=F, sep="\t", row.names=F)
}

message("completed")
