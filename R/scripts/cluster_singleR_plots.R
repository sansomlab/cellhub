## Run SingleR in a given Seurat object

# Libraries ----

stopifnot(
  require(optparse),
  require(dplyr),
  require(cellhub),
  require(SingleR),
  require(ggforce),
  require(MASS),
  require(viridis)
)

# Options ----

option_list <- list(
  make_option(c("--metadata"), default="metadata.tsv.gz",
              help="the cell metadata"),
  make_option(c("--scores"), default = NULL,
              help="tsv file with the singleR scores matrix"),
  make_option(c("--labels"), default = NULL,
              help="tsv file with the singleR labels"),
make_option(c("--reference"), default = "none",
              help="the name of the reference"),
  make_option(c("--show_annotation_in_plots"), default=NULL,
              help="Column names from the metadata slot to show in the annotation of the output plots"),
  make_option(c("--outdir"), default="seurat.out.dir",
              help="outdir"),
  make_option(c("--pdf"), default=FALSE,
              help="Should pdf plots be made?")

)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

message("reading metadata")
metadata <- read.table(opt$metadata, header=T, sep="\t",
                       row.names="barcode_id")

message("reading scores")
scores <- read.table(opt$scores, header=T, sep="\t",
                       row.names="barcode_id")

message("reading predictions]")
predictions <- read.table(opt$labels, header=T, sep="\t",
                          row.names="barcode_id")

scores <- scores[rownames(metadata), ]
predictions <- predictions[rownames(metadata), ]

scores$barcode_id <- NULL

predictions$scores <- as.matrix(scores)

# Plot label assignment scores
cat("Plotting prediction results \n")

qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mito")

if (!is.null(opt$show_annotation_in_plots)) {
    meta_cols <- unlist(strsplit(opt$show_annotation_in_plots, split = ","))
      keep_annotation <- c(qc_cols, meta_cols)
} else {
    keep_annotation <- qc_cols
}

    keep_annotation <- keep_annotation[keep_annotation %in% colnames(metadata)]

metadata <- metadata[,keep_annotation]

if(nrow(predictions) > 50000)
{
    message("downsampling to 50,000 cells")
    cells.use = sample(rownames(predictions), 50000)
} else {
    cells.use = rownames(predictions)
}

message("making the plot")
do_plot <- function() {
    plotScoreHeatmap(predictions, show.labels = TRUE,
                     cells.use=cells.use,
                     annotation_col=metadata)
}

save_plots(file.path(opt$outdir, 
                     paste(opt$reference, "heatmap", sep=".")),
           plot_fn=do_plot,
           width = 12,
           to_pdf=opt$pdf,
           height = 8)

message("Completed")
