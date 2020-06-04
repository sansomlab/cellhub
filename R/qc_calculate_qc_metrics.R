# Libraries ------
stopifnot(require(optparse), 
          require(futile.logger),
          require(dplyr),
          require(DropletUtils),
          require(tibble),
          require(scater),
          require(Matrix)
)

# Global options ------

options(stringsAsFactors = F)


# Script arguments ------

option_list <- list(
  make_option(c("--cellranger_dir"), default=".",
              help="Folder with cellranger output. Must include barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--genesets_file"), default=NULL,
              help="Two-column tsv file with genesets to evaluate, one geneset per row. First column: name of geneset; Second column: name of the file containing the geneset. The file containing the geneset must be a header-less, one-column tsv file with gene names (one per row)."),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to run scater's perCellQCMetrics function"),
  make_option(c("--log_filename"), default="qc_metrics.log"),
  make_option(c("--outfile"), default="qc_metrics.tsv.gz")
)

opt <- parse_args(OptionParser(option_list=option_list))


# Logger ------

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)


# Number of cores to use perCellQC function ------

multicoreParam <- MulticoreParam(workers = opt$numcores)


# Read in data ------
#####################

flog.info("Reading data and making SingleCellExperiment object...")
s <- read10xCounts(file.path(opt$cellranger_dir))
flog.info("Number of cells in input: %s", format(ncol(s), big.mark=","))


# Compute basic QC metrics ------
#################################

flog.info("Calculating basic QC metrics...")

# Number of genes
ngenes <- colSums(counts(s) > 0)

# Total UMI counts
total_UMI <- colSums(counts(s))

# Create dataframe 
cell_qc <- tibble(barcode = colData(s)$Barcode, 
                  ngenes = ngenes, 
                  total_UMI= total_UMI) 


# Compute percentage of genesets ------
#######################################

flog.info("Computing percentage of genesets...")

flog.info("Creating list of genesets to query")
genesets <- list()

# Mitochondrial genes
genesets$pct_mitochondrial <- grep("^MT-", rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s mitochondrial genes for pct_mitochondrial", 
          length(genesets$pct_mitochondrial))

# Ribosomal protein genes
genesets$pct_ribosomal <- grep("^RPS|^RPL", rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s ribosomal protein genes for pct_ribosomal", 
          length(genesets$pct_ribosomal))

# Immunoglobin genes
genesets$pct_immunoglobin <- grep("^IGH|^IGL|^IGK", rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s immunoglobin genes for pct_immunoglobin", 
          length(genesets$pct_immunoglobin))

# Hemoglobin genes
genesets$pct_hemoglobin <- grep("^Hbb|^Hba|^Hbq", rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s immunoglobin genes for pct_hemoglobin", 
          length(genesets$pct_hemoglobin))

# Add any other genesets provided
if (!is.null(opt$genesets_file)){
  
  flog.info("Reading genesets file %s", opt$genesets_file)
  provided_genesets <- read.csv(opt$genesets_file, sep="\t", header=FALSE)
  
  for (i in 1:nrow(provided_genesets)) {
    geneset_name <- provided_genesets[i,1]
    
    # Read geneset file
    flog.info("Reading %s geneset file", geneset_name)
    geneset <- read.csv(provided_genesets[i,2], sep="\t", header=FALSE)[,1]
    flog.info("Geneset %s contains %s genes", 
              geneset_name, format(length(geneset), big.mark=","))
    
    # Add to 'genesets' list
    genesets[[geneset_name]] <- which(rowData(s)$Symbol %in% geneset)
    flog.info("Including %s genes for %s", 
              format(length(geneset), big.mark=","), geneset_name)
  }
}

flog.info("Running scater function for computing percentage of genesets...")
pct_genesets <- perCellQCMetrics(s, subsets=genesets,
                                 BPPARAM = multicoreParam)

# Format output
pct_genesets <- pct_genesets[,grep("_percent", colnames(pct_genesets))]
colnames(pct_genesets) <- gsub("^subsets_", "", colnames(pct_genesets))
colnames(pct_genesets) <- gsub("_percent$", "", colnames(pct_genesets))
pct_genesets <- pct_genesets %>% as_tibble()
pct_genesets$barcode <- colData(s)$Barcode

# Append to QC table
cell_qc <- left_join(cell_qc, pct_genesets, by = "barcode")


# Compute mitoribo ratio ------
###############################

# create matrix with gene symbol as rownames
m <- counts(s)
rnames <- left_join(data.frame(ID=rownames(m)), 
                    as.data.frame(rowData(s)),
                    by="ID")
if (length(rownames(m)) != nrow(rnames)) {
  stop("Rownames in matrix do not have perfect matches to gene symbols")
}
rownames(m) <- rnames$Symbol

# Calculate mitoribo ratio
flog.info("Calculating mitoribo ratio...")
GetMitRibRatio <- function(matrix) {
  mitrib <- grep(pattern = "^MT-|^RP[SL]", x = rownames(x = matrix), value = TRUE)
  mito <- grep(pattern = "^MT-", x = rownames(x = matrix), value = TRUE)
  mitribcounts<- matrix[which(rownames(matrix) %in% mitrib), ]
  mitoribo_ratio <- Matrix::colSums(mitribcounts[mito, , drop = FALSE])/Matrix::colSums(mitribcounts)
  mitoribo_ratio <- as.data.frame(mitoribo_ratio)
  return(mitoribo_ratio)
}
mitoribo_ratio <- GetMitRibRatio(m)
mitoribo_ratio$barcode <- colData(s)$Barcode

# Append to QC table
cell_qc <- left_join(cell_qc, mitoribo_ratio, by = "barcode")

# Write table
flog.info("Writing output table")
write.table(cell_qc, file = gzfile(opt$outfile), quote=FALSE, row.names = FALSE)

flog.info("Completed")