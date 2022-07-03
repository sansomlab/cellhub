## Title ----
##
## Fetch geneset annotations.
##
## Description ----
##
## This script retrieves ID mappings (Ensembl to entrez).
##
## Details ----
##
## Ensembl to Entrez mappings are retrieved using biomaRt.
##
## Usage ----
##
## $ Rscript getGenesetAnnotations.R
##           --ensemblversion=latest
##           --species=mm
##           --outdir=.


# Libraries ----

stopifnot(
  require(optparse),
  require(gsfisher)
)

# Options ----

option_list <- list(
    make_option(
      c("--ensemblversion"),
      default="latest",
      help="either latest or a specific number"
    ),
    make_option(
      c("--ensemblhost"),
      default=NULL,
      help="the ensembl host address"
    ),
    make_option(
      c("--species"),
      default="none",
      help="species - mm or hs"
      ),
    make_option(
      c("--outfile"),
      default="none",
      help="outfile")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Fetch Ensembl to Entrez ID mappings ----
if(opt$ensemblversion=="latest")
{
    version <- NULL
} else { version <- opt$ensemblversion }

if(is.null(opt$ensembl_host))
{
    # use the default host (www.ensembl.org)
    anno <- fetchAnnotation(species=opt$species,
                            ensembl_version=version)
} else {
    anno <- fetchAnnotation(species=opt$species,
                            ensembl_version=version,
                            ensembl_host=opt$ensemblhost)
}

write.table(anno,
            gzfile(opt$outfile),
            quote=FALSE,
            row.names=FALSE,sep="\t")

