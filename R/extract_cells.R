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
      c("--outdir"),
      default=".",
      help="The name of the directory for the output"
      )
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
message(opt)

## Read in the matrix
matrix_path <- file.path(opt$matrixdir, "matrix.mtx.gz")
matrix_barcodes <- file.path(opt$matrixdir, "barcodes.tsv.gz")
matrix_features <- file.path(opt$matrixdir, "features.tsv.gz")

x <- readMM(matrix_path)
colnames(x)  <- read.table(matrix_barcodes, stringsAsFactors=FALSE)$V1
features <- read.table(matrix_features, stringsAsFactors=FALSE, sep = "\t")
rownames(x)  <- features$V1

message("matrix dimensions before subsetting")
message(dim(x))

## Read in the list of cells to extract
cells <- read.table(gzfile(opt$cells),stringsAsFactors=FALSE)$V1
message("\n number of cells to extract \n")
message(length(cells))
message("\nnumber of cells to extract present in ref mtx \n")
message(length(which( colnames(x)%in%cells)))

x <- x[, colnames(x)%in%cells]
out_folder = opt$outdir
message("matrix dimensions after subsetting")
message(dim(x))

## Save the results
matrix_path = file.path(out_folder, "matrix.mtx")
barcodes_path = file.path(out_folder, "barcodes.tsv")
features_path = file.path(out_folder, "features.tsv.gz")

## write out the matrix
writeMM(x, matrix_path)

## write out the "cell" barcodes
write.table(colnames(x),
    	    barcodes_path,
            col.names=FALSE,
            sep=",",
            row.names=FALSE,
            quote=FALSE)

## copy over the features file
file.copy(matrix_features, features_path)

## compress the outputs
gzip(barcodes_path, overwrite=TRUE)
gzip(matrix_path, overwrite=TRUE)
