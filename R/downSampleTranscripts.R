## Title ----
##
## Subset the aggregated gene exp. matrix to equalize seq-depth differences among 
## given groups of cells
##
## Description ----
##
## Subset the aggregated gene exp. matrix to equalize seq-depth differences among 
## given groups of cells. It accepts the outputs of the pipeline_fetch_cells.py
##
## The script can downsample UMI counts across a groups of cells
## to match their median or mean.
##
## Usage ----
##
## statement = '''Rscript %(tenx_dir)s/R/downSampleTranscripts.R
##                --cellgroupvar=<string> (e.g. sample_id, batch_id, etc.)
##                --genexpdir=<path> (e.g. output.dir)
##                --downsample=<string> (e.g. median)
##                --outdir=<path> (e.g. %(out_dir)s )
##                &> %(outfile)s
##             '''

message("downSampleTranscripts.R here we go!")
timestamp()

# Libraries ----

stopifnot(
  require(optparse),
  require(methods), # https://github.com/tudo-r/BatchJobs/issues/27
  require(Matrix),
  require(S4Vectors),
  require(ggplot2),
  require(R.utils),
  require(DropletUtils),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
        c("--genexpdir"),
        help="matrix directory"
    ),
    make_option(
        c("--cellmetadata"),
        help="Path to tsv.gz file with cell-metadata"
    ),
    make_option(
        c("--cellgroupvar"),
        help="cell-grouping variable to equalize seq-depth"
    ),
    make_option(
        c("--downsample"),
        default="median",
        help="Strategy to normalise UMI between cell groups"
    ),
    make_option(
        c("--outdir"),
        default=".",
        help="location where the outputs will be written"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Functions ----

downsampleMatrix <- function(matrixUMI, downsample_method="median", library_ids=c())
{
  # downsampleCounts drops the colnames
  backup_colnames <- colnames(matrixUMI)

  # Total UMI count per cell
  nUMIs <- Matrix::colSums(matrixUMI)

  # Summary UMI metric by sample (e.g. median)

  statUMIs <- tapply(nUMIs, library_ids, downsample_method)
  cat(sprintf("%s(UMIs) before downsampling (by agg_id):\n", downsample_method))
  print(statUMIs)

  # Identify proportion to retain for each sample
  normFactorUMIs <- min(statUMIs) / statUMIs
  cat("Proportion for downsampling by agg_id:\n")
  print(normFactorUMIs)

  # Extend scaling factor to all cells by sample
  colNormFactor <- normFactorUMIs[library_ids]

  # Downsample
  matrixUMI <- DropletUtils::downsampleMatrix(as(matrixUMI, "dgCMatrix"), colNormFactor, bycol=TRUE)

  # downsampleCounts drops the colnames
  colnames(matrixUMI) <- backup_colnames

  # Ensure the downsampling has worked properly
  nUMIs <- Matrix::colSums(matrixUMI)
  statUMIs <- tapply(nUMIs, library_ids, downsample_method)
  cat(sprintf("%s(UMIs) after downsampling (by agg_id):\n", downsample_method))
  print(statUMIs)
  matrixUMI
}

plotDownsamp <- function(matrixUMI, metadata, basename) {

    # Total UMI count per cell
    nUMIs <- Matrix::colSums(matrixUMI)

    # Collate total UMI and cell rank for each cell of each sample
    inputStats <- do.call("rbind", lapply(
        unique(metadata$sample_id),
        function(id){
            codes_id <- subset(metadata, sample_id == id, "barcode", drop=TRUE)
            nUMIs_id <- nUMIs[codes_id]
            data.frame(
                "nUMIs"=sort(nUMIs_id, decreasing=TRUE),
                "CellRank"=seq_along(nUMIs_id),
                "sample"=id
            )
        }))

    message(head(inputStats, 3))
    print(head(inputStats, 3))

    # UMI vs Rank ----

    maxCellRank <- max(inputStats$CellRank)

    gg <- ggplot(inputStats) +
        geom_line(aes(x = `CellRank`, y = nUMIs, colour=sample), size=0.25) +
        scale_x_log10(
            limits=c(1, maxCellRank)
        ) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs)))) 
        ) +
        annotation_logticks() +
        labs(y="UMI count", x="Cell rank") +
        theme_bw() +
        theme(
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank()
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_ranked.pdf", basename)), gg,
        width=7, height=5
    )

    # UMI violion plot ----

    gg <- ggplot(inputStats) +
        geom_violin(
            aes(x = sample, y = nUMIs, colour=sample),
            draw_quantiles=c(0.5), size=0.25) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs)))) 
        ) +
        annotation_logticks(sides="l") +
        labs(y="UMI count", x="Sample") +
        theme_bw() +
        theme(
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(angle=90)
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_violin.pdf", basename)), gg,
        width=7, height=5
    )

    return(TRUE)
}

downsample <- function(matrixUMI, metadata, agg_ids)
{

  message("making pre-downsampling diagnostic plot")
  plotDownsamp(matrixUMI, metadata, "UMI_input")

  message("downsampling the matrix")
  matrixUMI <- downsampleMatrix(matrixUMI,
                                downsample_method=opt$downsample,
                                library_ids=agg_ids)

  message("making post-downsampling diagnostic plot")
  plotDownsamp(matrixUMI, metadata, "UMI_downsampled")
  matrixUMI
}


# Setup output directory ----

if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir)
}

# Input data ----

## Matrix

matrixFile <- file.path(opt$genexpdir, "matrix.mtx.gz")
stopifnot(file.exists(matrixFile))
cat("Importing matrix from:", matrixFile, " ... ")
matrixUMI <- readMM(gzfile(matrixFile))

cat("Done.\n")
cat(
    "Input matrix size:",
    sprintf("%i rows/genes, %i columns/cells\n", nrow(matrixUMI), ncol(matrixUMI))
)

## Barcodes

barcodeFile <- file.path(opt$genexpdir, "barcodes.tsv.gz")
stopifnot(file.exists(barcodeFile))
cat("Importing cell barcodes from:", barcodeFile, " ... ")
barcodes <- scan(gzfile(barcodeFile), "character")

## Metadata

metadataFile <- file.path(opt$cellmetadata)
stopifnot(file.exists(metadataFile))

cat("Importing metadata from:", metadataFile, " ... ")
metadata <- data.frame(data.table::fread(metadataFile))
cat("Done.\n")
cat(
    "Input metadata size:",
    sprintf("%i rows/cells, %i columns\n", nrow(metadata), ncol(metadata)),
    "\n"
)

# Preprocess ----

## Metadata
metadata$agg_id <- as.factor(metadata[[opt$cellgroupvar]])
metadata$barcode <- as.character(metadata$barcode_id)

## UMI matrix colnames / cell barcodes
colnames(matrixUMI) <- barcodes

cat("Applying downsampling with no-subsetting.\n")

matrixUMI <- downsample(matrixUMI, metadata, metadata$agg_id)

# Print information
cat(
    "New matrix size (cells subset):",
    sprintf("%i rows/genes, %i columns/cells\n", nrow(matrixUMI), ncol(matrixUMI))
)

# Write out matrix ----

genesFile <- file.path(opt$genexpdir,"features.tsv.gz")
writeMatrix(opt$outdir, matrixUMI, barcodes, genesFile, metadata)

timestamp()
message("Completed")
