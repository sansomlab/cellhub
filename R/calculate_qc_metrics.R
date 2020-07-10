# Libraries ------
stopifnot(require(optparse), 
          require(futile.logger),
          require(dplyr),
          require(DropletUtils),
          require(tibble),
          require(scater),
          require(Matrix),
          require(cellqc)
)

# Global options ------

options(stringsAsFactors = F)


# Script arguments ------

option_list <- list(
  make_option(c("--cellranger_dir"), default=".",
              help="Folder with filtered cellranger output. Must include barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--sample"), default=NULL,
              help="sample or channel id"),
  make_option(c("--barcodes_to_label_as_True"), default=NULL,
              help="A header-less, two-column, gziped tsv.gz file with barcodes to label. 
              The name of the file will be used to set the column name in the output dataframe ('.tsv' will be removed from the final name).
              First column: barcode name in the format UMI-1, e.g AAACCTGAGAGGTACC-1.
              Second column: sample or channel id name.
              If a barcode is included in the input file, it will be labeled as 'True', otherwise as 'False'"),
  make_option(c("--genesets_file"), default=NULL,
              help="Two-column tsv file with genesets to evaluate, one geneset per row. 
              First column: name of geneset; 
              Second column: name of the file containing the geneset. 
              The file containing the geneset must be a header-less, one-column tsv file with gene names (one per row)."),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to run scater's perCellQCMetrics function"),
  make_option(c("--log_filename"), default="qcmetrics.log"),
  make_option(c("--outfile"), default="qcmetrics.tsv.gz")
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
                  sample=opt$sample,
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
genesets$pct_immunoglobin <- grep("^IGH|^IGLC|^IGLJ|^IGLL|^IGLV|^IGK", 
                                  rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s immunoglobin genes for pct_immunoglobin", 
          length(genesets$pct_immunoglobin))

# Hemoglobin genes
genesets$pct_hemoglobin <- grep("^Hbb|^Hba|^Hbq", rowData(s)$Symbol, ignore.case=TRUE)
flog.info("Including %s hemoglobin genes for pct_hemoglobin", 
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
              format(length(genesets[[geneset_name]]), big.mark=","), geneset_name)
    
    # Report genes not found
    if ( sum(!geneset %in% rowData(s)$Symbol) >= 1) {
      not_found <- geneset[-which(geneset %in% rowData(s)$Symbol)]
      not_found <- paste0(not_found, collapse = ", ")
      flog.info("Genes not found in count matrix: %s", 
                not_found)
    }
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
mitoribo_ratio <- GetMitRibRatio(m)
mitoribo_ratio$barcode <- colData(s)$Barcode

# Append to QC table
cell_qc <- left_join(cell_qc, mitoribo_ratio, by = "barcode")


# Label barcodes ------
#######################
if (!is.null(opt$barcodes_to_label_as_True)){
  bc <- read.csv(gzfile(opt$barcodes_to_label_as_True), sep="\t", header=FALSE) 
  colnames(bc) <- c("barcode", "sample")
  # Filter barcodes from current sample
  bc <- bc[bc$sample %in% opt$sample,]
  # Label barcodes in output dataframe
  cell_qc$newcol <- "False"
  cell_qc$newcol[cell_qc$barcode %in% bc$barcode] <- "True"
  # Set column name in output
  colname <- opt$barcodes_to_label_as_True
  colname <- gsub(".tsv.gz", "", gsub(".*/", "", colname))
  colnames(cell_qc)[colnames(cell_qc) == "newcol"] <- colname
}

# Format ouput to match tables in other pipelines from the cellhub repository
colnames(cell_qc)[colnames(cell_qc)=="barcode"] <- "BARCODE"

# Write table
flog.info("Writing output table")
write.table(cell_qc, file = gzfile(opt$outfile), quote=FALSE, row.names = FALSE, sep="\t")

flog.info("Completed")