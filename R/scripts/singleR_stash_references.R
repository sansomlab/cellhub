# Libraries ----

stopifnot(
  require(celldex),
  require(ensembldb),
  require(optparse)
)

# Options ----

option_list <- list(
  make_option(c("--outdir"), default=".",
              help="outdir")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:\n")
print(opt)

references <- c("HumanPrimaryCellAtlasData",
                "BlueprintEncodeData",
                "MouseRNAseqData",
                "ImmGenData",
                "DatabaseImmuneCellExpressionData",
                "NovershternHematopoieticData",
                "MonacoImmuneData")

for(ref in references)
{
  message("Fetching ", ref)
  ref.x <- match.fun(ref)
  ref.x <- ref.x(ensembl=TRUE)
  
  message("Stashing ", ref)
  saveRDS(ref.x, 
          file.path(opt$outdir, paste(ref, "rds", sep=".")))
  
}

message("Completed")
