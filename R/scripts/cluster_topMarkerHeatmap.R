## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(openxlsx),
  require(data.table),
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
    make_option(c("--clusterids"), default="none",
                help="A tsv file containing the cluster identities"),
    make_option(c("--metadata"), default="none",
                help="A tsv file containing the metadata"),
    make_option(c("--markers"), default="none",
                help="The table of summarised markers"),
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


message("Making a heatmap of the top marker genes from each cluster")

markers <- read.table(opt$markers,
                      sep="\t", header=T, as.is=T)

filtered_markers <- data.table(markers[markers$p.adj < 0.1 & !is.na(markers$p.adj),])

## make a heatmap of the top DE genes.
# filtered_markers %>% group_by(cluster) %>% top_n(20, log2FC) -> top20

# if(!is.null(opt$subgroup))
# {
#     if(!opt$subgroup %in% colnames(s@meta.data))
#     {
#         opt$subgroup <- NULL
#     }
# }

if(!is.null(opt$subgroup)) { opt$subgroup <- strsplit(opt$subgroup,",")[[1]]}

mch <- markerComplexHeatmap(loom_path=opt$loom,
                            matrix_loc=opt$matrix_loc,
                            barcode_id_loc=opt$barcode_id_loc, #"col_attrs/barcode_id",
                            gene_id_loc=opt$gene_id_loc, # "row_attrs/gene_ids",
                            scale=opt$scale,
                            cluster_ids=opt$clusterids,
                            metadata_file=opt$metadata,
                            marker_table=filtered_markers,
                            priority="pval",
                            n_markers=20,
                            cells_use=NULL,
                            row_names_gp=10,
                            sub_group=opt$subgroup)

drawHeatmap <- function()
{
    draw(mch)
}

save_plots(paste(outPrefix,"heatmap", sep="."),
           plot_fn=drawHeatmap,
           width = 7,
           to_pdf = opt$pdf,
           height = 9)

message("completed")
