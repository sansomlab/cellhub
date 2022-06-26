## make seurat object from anndata output

stopifnot(require(Seurat))
stopifnot(require(yaml))
stopifnot(require(rhdf5))
stopifnot(require(optparse))
stopifnot(require(futile.logger))
stopifnot(require(SeuratDisk))

# these parameters are passed from Rscript run
option_list <- list(
  make_option(
    c("--task_yml"),
    dest = "task_yml",
    help="Path to yml file"
  ),
  make_option(
    c("--log_filename"),
    dest = "log_filename",
    help="Path to log file"
  ))
params <- parse_args(OptionParser(option_list=option_list))

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))
flog.info("Running with parameters: %s", params, capture = TRUE)


# these options are read from the yml
default_options <- list(
  # The path to the directory containing the
  # matrix and associates barcodes and features
  "matrixdir" = NULL,

  # Path to the output folder
  "outdir" = "sample.dir",

  # Which format to write? Options rds|h5seurat
  "output_format" = "rds",


  # Only keep seurat data slot? NULL = keep all, 
  # 'reduce' = keep only data
  "seurat_data_only" = NULL,
 
  # Name for the new dim reduction in seurat object
  # and file name within the folder for infile_h5
  "dim_name" = "harmony",

  # cell numbers for subsets
  "cell_numbers" = NULL

)

options <- read_yaml(params$task_yml)

# Update the default options
if(!is.null(options)) {
  opt <- utils::modifyList(default_options, options)
} else{
  opt <- default_options
}

flog.info("Running with options: ", opt, capture = TRUE)
flog.info("\n")

if (is.null(opt$matrixdir)){
   h5ad_input = file.path(opt$outdir, "full_anndata.h5ad")
   h5seurat_input = gsub(".h5ad", ".h5seurat", h5ad_input)
   flog.info("Convert anndata to h5seurat")
   Convert(h5ad_input, dest = "h5seurat", overwrite = TRUE)

   # load anndata with raw counts (=counts) and log-norm slot (=data)
   s <- LoadH5Seurat(h5seurat_input)
} else {
  data <- Read10X(opt$matrixdir)
  s <- CreateSeuratObject(counts=data,
                          project="COMBAT_50k_subsample"
			  )
  s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
  metafile <- file.path(opt$outdir, "metadata.tsv.gz")

  metadata <- read.table(gzfile(metafile),
  	                        sep="\t", header=TRUE, as.is=TRUE)
  rownames(metadata) <- metadata$barcode
  metadata$barcode <- NULL
  metadata <- metadata[colnames(x = s), ]
  for(meta_col in colnames(metadata)) {
  s[[meta_col]] <- metadata[[meta_col]]
  }

}

flog.info("Read dim reduction")

h5file = H5Fopen(file.path(opt$outdir, "scaled.h5"))

#dim_reduction, dim_name
#comp <- read.table(gzfile(file.path(folder_run, opt$dim_reduction)), sep = "\t", header=TRUE)
comp <- t(h5file$comp)
rownames(comp) = h5file$barcodes
colnames(comp) = paste("harmony", c(1:ncol(comp)), sep="_")

s[[opt$dim_name]] <- CreateDimReducObject(embeddings = as.matrix(comp), assay = "RNA")

if (!is.null(opt$seurat_data_only)){
  flog.info("Do not add scaled data, but instead run DietSeurat before writing out the object (below)")
} else {
  flog.info("Read in scaled data")

  scaled_data = h5file$scaled_data
  colnames(scaled_data) = h5file$barcodes
  rownames(scaled_data) = h5file$genes

  scaled_data = scaled_data[rownames(scaled_data) %in% rownames(s),]
  # set scaled data assay
  s <- SetAssayData(object = s, slot = "scale.data",
      new.data = scaled_data, assay = "RNA")
}

VariableFeatures(s) <- h5file$genes

## change all factors to characters as otherwise
## not converted to anndata properly!

flog.info("Fixing factors in metadata")
metadata_new = s[[]]
i <- sapply(metadata_new, is.factor)
metadata_new[i] <- lapply(metadata_new[i], as.character)
s <- AddMetaData(s, metadata=metadata_new)

# the gene_ids are also factors
if (ncol(s@assays$RNA@meta.features) > 0){
   flog.info("Fixing factors in meta.features slot")
   meta_features = s@assays$RNA@meta.features
   i <- sapply(meta_features, is.factor)
   meta_features[i] <- lapply(meta_features[i], as.character)
   s@assays$RNA@meta.features <- meta_features
}


if(!is.null(opt$cell_numbers)){
    flog.info("Start making subsets")
    for (n in opt$cell_numbers){
     subset_seurat =  subset(s, cells = sample(Cells(s), size=as.numeric(n)))
     flog.info("starting to write seurat object")
     if (opt$output_format == "rds"){
       saveRDS(subset_seurat, file=file.path(opt$outdir, paste0("seurat_object_", n, ".rds")))
     } else {
       SaveH5Seurat(subset_seurat, file=file.path(opt$outdir,
       				      paste0("seurat_object_", n, ".h5seurat")),
				      overwrite = TRUE)
     }
     flog.info("finished writing object for %s cells", dim(subset_seurat)[2])
    }
}

if (!is.null(opt$seurat_data_only)){
   flog.info("Run DietSeurat")
   s <- DietSeurat(s, counts=FALSE, data=TRUE, scale.data = FALSE, dimreducs= opt$dim_name)
}

flog.info("Start writing main seurat object")

if (opt$output_format == "rds"){
    saveRDS(s, file.path(opt$outdir, "seurat_object.rds"))
} else {
    SaveH5Seurat(s, file.path(opt$outdir, "seurat_object.h5seurat"),
    		 overwrite = TRUE)
}

flog.info("Completed")
