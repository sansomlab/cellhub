# Libraries ----

stopifnot(
  require(SingleR),
  require(DropletUtils),
  require(celldex),
  require(ensembldb),
  require(BiocParallel),
  require(optparse)
)

# Options ----

option_list <- list(
  make_option(c("--mtxdir"), default="matrix.dir",
              help="The location of the 10x mtx directory"),
  make_option(c("--sample"), default="none",
              help="The name of the sample"),
  make_option(c("--refstashdir"),
              help="folder with the rds stash of reference datasets" ),
    make_option(c("--reference"),
              help="the name of the celldex reference" ),
  make_option(c("--workers"), default=10,
              help="number of parallel processes to use"),
  make_option(c("--outdir"), default="seurat.out.dir",
              help="outdir")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:\n")
print(opt)

message("Reading in the count data")
sce <- read10xCounts(opt$mtxdir)

# sample_suffix = paste0("-", opt$sample)
# colnames(sce) <- gsub("-1$", sample_suffix, colData(sce)$Barcode)
colnames(sce) <- colData(sce)$Barcode

message("Loading the reference dataset: ", opt$reference)

ref_path = file.path(opt$refstashdir, paste(opt$reference,"rds",sep="."))
print(ref_path)
ref.se <- readRDS(ref_path)

if (grepl("ENSMUSG", rownames(ref.se)[1]) & grepl("ENSG",rownames(sce)[1]))
{
  message("cleaning up mouse ensembl indentifiers")
  rownames(ref.se) <- gsub("ENSMUSG", "ENSG", rownames(ref.se))
}
message("setting up the parallel environment")
multicoreParam <- MulticoreParam(workers = opt$workers)

message("running singleR")
pred <- SingleR(test = sce,
                ref = ref.se,
                assay.type.test = "counts",
                labels = ref.se$label.main,
                BPPARAM = multicoreParam)

message("saving the labels")
labels <- data.frame(barcode=rownames(pred),
                     library_id=rep(opt$sample,nrow(pred)),
                     labels=pred$labels,
                     pruned.labels=pred$pruned.labels)

write.table(labels,
            gzfile(file.path(opt$outdir,paste(opt$sample,"labels.tsv.gz", sep="."))),
            col.names=TRUE, row.names=FALSE,
            sep="\t", quote=FALSE)

message("saving the scores")
scores <- data.frame(barcode=rownames(pred),
                     library_id=rep(opt$sample,nrow(pred)))
scores <- cbind(scores, data.frame(pred$scores))

write.table(scores,
            gzfile(file.path(opt$outdir,paste(opt$sample,"scores.tsv.gz", sep="."))),
            col.names=TRUE, row.names=FALSE,
            sep="\t", quote=FALSE)

message("Completed\n\n")
