## Title ----
##
## Fetch geneset annotations.
##
## Description ----
##
## This script retrieves KEGG pathway information.
##
## Details ----
##
## KEGG pathways are retrieved
## directly from KEGG.
##
## Usage ----
##


# Libraries ----

stopifnot(
  require(optparse),
  require(gsfisher)
)

# Options ----

option_list <- list(
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

# Fetch KEGG pathways ----
kegg_pathways <- fetchKEGG(species=opt$species)
saveRDS(kegg_pathways, file=file.path(opt$outfile))
