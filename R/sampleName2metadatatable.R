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
#' 
#' This script will get the output of pipeline_cellranger_multi.py in the 
#' folder ./data.dir to assembly a input_samples.tsv file containing the
#' available sample sequencing metadata provided in the configuration *.yml
#' file of the pipeline_cellranger_multi.py
#' ---

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

samples_df <- lapply(samples, function(s) {
  
  sample_id <- gsub(".csv", "", basename(s))
  smeta <- readLines(paste0("data.dir/", s))
  smeta <- smeta[sapply(smeta, function(x) x != "")]
  sample_dat <- lapply(seq_along(smeta), function(i) {
    
    sample_meta <- smeta[i]
    if(!grepl("\\[", sample_meta)) {
      
      infos <- unlist(strsplit(sample_meta, ","))
      
      if(length(infos) == 2) {
        
        df <- data.frame('i' = infos[2])
        colnames(df) <- infos[1]
        df
        
      } else if(length(infos) > 2 & i < length(smeta)) {
        
        cols <- infos
        cont <- unlist(strsplit(smeta[i + 1], ','))
        cont <- t(data.frame(cont))
        colnames(cont) <- cols[1:ncol(cont)]
        cont
        
      }
    }
  })
  
  sample_dat <- Filter(Negate(is.null), sample_dat)
  sample_df <- do.call(cbind, sample_dat)
  sample_df[["sample_id"]] <- sample_id
  sample_df[["filt_path"]] <- paste0(getwd(), '/',
                                     sample_id,
                                     "/outs/per_sample_outs/",
                                     sample_id,
                                     "/count/sample_feature_bc_matrix")
  sample_df[["raw_path"]] <- paste0(getwd(), '/',
                                    sample_id,
                                    "/outs//multi/count/raw_feature_bc_matrix")
  sample_df
  
})

input_sample_meta <- do.call(rbind, samples_df)

write.table(input_sample_meta, params$outfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat(paste0("Sample metadata ", params$outfile), "has been created.")