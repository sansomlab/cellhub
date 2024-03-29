## Cluster the single cells in a given Seurat object

# Libraries ----

stopifnot(
  require(optparse),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(cellhub),
  require(xtable)
)

# Options ----

option_list <- list(

    make_option(c("--clusters"), default="scanpy.clusters.tsv.gz",
                help="the scanpy cluster assignments"),
    make_option(c("--predefined"), default="none",
                help="a file containing a set of predefined clusters"),
    make_option(c("--mincells"), type="integer", default=10,
                help="clusters with fewer cells are set to NA"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


if(opt$predefined=="none")
    {
        clusters <- read.table(opt$clusters, sep="\t", header=T, as.is=T)
        cluster_ids <- clusters$cluster_id
        x <- table(cluster_ids)

        rejected_clusters <- names(x[x<opt$mincells])
        cluster_ids[cluster_ids %in% rejected_clusters] <- "911"


    } else {

        clusters <- read.table(opt$predefined, sep="\t", header=T, as.is=T)
        cluster_ids <- clusters$cluster_id

        }

cluster_ids <- factor(as.numeric(cluster_ids))
names(cluster_ids) <- clusters$barcode_id

unique_cluster_ids <- unique(cluster_ids)
print(unique_cluster_ids)

message(sprintf("saving a unique list of cluster ids"))
write.table(unique_cluster_ids, file=file.path(opt$outdir,"cluster_ids.tsv"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)



cluster_colors <- gg_color_hue(length(unique_cluster_ids))
message(sprintf("saving the cluster colors (ggplot)"))
write.table(cluster_colors, file=file.path(opt$outdir,"cluster_colors.tsv"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)

cluster_assignments <- data.frame(barcode_id=as.character(names(cluster_ids)),
                                  cluster_id=as.character(cluster_ids))

write.table(cluster_assignments,
            gzfile(file.path(opt$outdir, "cluster_ids.tsv.gz")),
            sep="\t", col.names=T, row.names=F, quote=F)



message(sprintf("saving numbers of cells per cluster"))
ncells_per_cluster <- table(cluster_assignments$cluster_id)

write.table(data.frame(cluster=names(ncells_per_cluster),
                       ncells=as.numeric(ncells_per_cluster)), 
            sep="\t",
            file=file.path(opt$outdir,"cluster_cell_counts.tsv"),
            quote=FALSE, col.names = TRUE, row.names = FALSE)

message("Completed")
