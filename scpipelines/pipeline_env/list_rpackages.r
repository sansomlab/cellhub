# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
}))

# Parse options ----

option_list <- list(
  make_option(
    c("--outfile", "-o"),
    default = 'env.out',
    help="Name of output file to append SessionInfo R packages")
)

opt <- parse_args(OptionParser(option_list = option_list))

message("Running with options:")
print(opt)

# Process data ----

sink(opt$outfile, append = TRUE)
cat("\n---R packages with sessioninfo---\n")
print(sessioninfo::session_info())
sink()

# Conclusion ---

message("Completed")
timestamp()
