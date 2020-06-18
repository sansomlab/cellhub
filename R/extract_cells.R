stopifnot(
  require(optparse),
  require(Matrix),
  require(R.utils)
)

# Options ----

option_list <- list(
    make_option(
      c("--cells"),
      default="none",
      help="a list of the cell barcodes to retain"
    ),
    make_option(
        c("--matrixdir"),
        default="none",
        help="the folder containing the input matrices"),
    make_option(
        c("--matrixid"),
        default="none",
        help="the name of subfolder containing the input matrix"),
    make_option(
      c("--matrixtype"),
      default="mm",
      help="the type of the input matrix, mm or loom"
      ),
    make_option(
      c("--outdir"),
      default=".",
      help="The name of the directory for the output"
      )
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

## Read in the matrix

outs = "outs/filtered_feature_bc_matrix/"

if(opt$matrixtype=="mm")
{
    matrix_path <- file.path(opt$matrixdir, opt$matrixid, outs, "matrix.mtx.gz")
    matrix_barcodes <- file.path(opt$matrixdir, opt$matrixid, outs, "barcodes.tsv.gz")
    matrix_features <- file.path(opt$matrixdir, opt$matrixid, outs, "features.tsv.gz")

    x <- readMM(matrix_path)
    colnames(x)  <- read.table(matrix_barcodes, stringsAsFactors=FALSE)$V1
} else if (opt$matrixtype=="loom") {

    stop("loom support not yet implemented")
} else {

    stop("matrix type not recognised")
}

print("matrix dimensions before subsetting")
print(dim(x))

## Read in the list of cells to extract
cells <- read.table(gzfile(opt$cells),stringsAsFactors=FALSE)$V1

print("number of cells to extract")
print(length(cells))

x <- x[, cells]

print("matrix dimensions after subsetting")
print(dim(x))


out_folder = opt$outdir

## Save the results
matrix_path = file.path(out_folder, "matrix.mtx")
barcodes_path = file.path(out_folder, "barcodes.tsv")
features_path = file.path(out_folder, "features.tsv.gz")

## write out the matrix
writeMM(x, matrix_path)

## write out the "cell" barcodes
write.table(read.table(text=colnames(x),sep="-")$V1,
  	    barcodes_path,
            col.names=FALSE,
            sep=",",
            row.names=FALSE,
            quote=FALSE)

## copy over the features file
file.copy(matrix_features, features_path)

## compress the outputs
gzip(matrix_path, overwrite=TRUE)
gzip(barcodes_path, overwrite=TRUE)
