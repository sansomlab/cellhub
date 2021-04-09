#' ---
#' title: "Sample names to sample metadata table"
#' output: 
#'  tsv_document
#' params:
#'  outfile:
#'  fig_path: 
#'  log_filename: "sample_design.metada.log"
#' ---
#' ---
#' Run makeSampleTable
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---

#+ setup, include=FALSE, echo=FALSE

timestamp()
message("Started")

# Libraries --------------------------------------------------------------------

stopifnot(require(dplyr),
          require(reshape2),
          require(optparse))

# Parameters -------------------------------------------------------------------

# these parameters are passed from Rscript run
option_list <- list(
  make_option(
    c("--outfile"),
    dest = "outfile",
    help="Path to output metadata file"
  ),
  make_option(
    c("--samplefiles"),
    dest = "samplenames",
    help = "list of input sample names"
  ))

params <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(params)

samples <- unlist(strsplit(params$samplenames, "///"))

outtab <- data.frame(
		     do.call(rbind, lapply(samples, function(s) {
					unlist(strsplit(s, "\\."))[1:3]
				     })
		             )
		     )

colnames(outtab) <- c("sample_id", "ncells", "exp_batch")

outtab[["filt_path"]] <- paste0(outtab[["sample_id"]], "-count/outs/filtered_feature_bc_matrix")
outtab[["raw_path"]] <- paste0(outtab[["sample_id"]], "-count/outs/raw_feature_bc_matrix")
outtab[["outs_path"]] <- paste0(outtab[["sample_id"]], "-count/outs")
outtab[["channel_id"]] <- 1:nrow(outtab)
outtab[["seq_batch"]] <- 1

write.table(outtab, params$outfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat(paste0("Sample metadata ", params$outfile), "has been created.")

