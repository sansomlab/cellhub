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
#' folder ./data.dir to assembly a input_libraries.tsv file containing the
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
      c("--libraryfiles"),
      dest = "librarynames",
      help = "list of input library names"),
  make_option(
      c("--librarydir"),
      default=".",
      dest = "librarydir",
      help = "list of input library names"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(opt)

libraries <- unlist(strsplit(opt$librarynames, ","))

libraries_df <- lapply(libraries, function(s) {

  library_id <- gsub(".csv", "", basename(s))
  smeta <- readLines(file.path(opt$librarydir, s))
  smeta <- smeta[sapply(smeta, function(x) x != "")]
  library_dat <- lapply(seq_along(smeta), function(i) {

    library_meta <- smeta[i]
    if(!grepl("\\[", library_meta)) {

      infos <- unlist(strsplit(library_meta, ","))

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

    library_dat <- Filter(Negate(is.null), library_dat)
    library_df <- do.call(cbind, library_dat)
    library_df[["library_id"]] <- library_id

    library_df[["filt_path"]] <- file.path(getwd(),
                                          opt$librarydir,
                                          library_id,
                                          "outs/per_sample_outs",
                                          library_id,
                                          "count/sample_feature_bc_matrix")
    library_df[["raw_path"]] <- file.path(getwd(),
                                         opt$librarydir,
                                         library_id,
                                         "outs/multi/count/raw_feature_bc_matrix")
  library_df

})

input_library_meta <- do.call(rbind, libraries_df)

write.table(input_library_meta,
            opt$outfile,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Library metadata ", opt$outfile, " has been created.")
