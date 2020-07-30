#' ---
#' title: "Ambient RNA summary"
#' output: 
#'  html_document:
#'   self_contained: false
#'   toc: true
#'   toc_float: true
#' params:
#'  task_yml: "/gfs/devel/mjgomez/tools/dropflow/Rmd/ambientRNA_summary.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "ambientRNA_summary.log"
#' ---
#' ---
#' Run ambient RNA summary
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(yaml),
          require(ggplot2),
          require(futile.logger),
          require(knitr),
          require(dplyr),
          require(ggrepel),
          require(stringr),
          require(tidyverse),
          require(ComplexHeatmap),
          require(RColorBrewer),
          require(grDevices),
          require(optparse))


# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE)

# Parameters -------------------------------------------------------------------

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

# The script expects the following paramters:
default_options <- list(
  # Path to the ambient RNA summary output file
  "outdir" = ".",
  
  # Comma separated paths to the directories containing the samples'
  # top genes data frame (top_ambient_genes.txt.gz), output of the ambient_rna.R
  "sample_indir" = "",
  
  # Comma separated names of the samples analysed, 
  # in the same order as the sample_indir
  "sample_id" = "",
  
  # File with sample metadata
  "sample_table" = ""
)


# here the yaml can also be read in from the default location in code 
# directory
options <- read_yaml(params$task_yml)

# Update the default options
if(!is.null(options)) {
  opt <- utils::modifyList(default_options, options)
} else{
  opt <- default_options
}


# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))

flog.info("Running with options: ", opt, capture = TRUE)
flog.info("\n")


## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

flog.info("Reading data of top ambient genes in samples")
indir <- str_trim(unlist(str_split(opt$sample_indir, pattern = ",")), 
                  side="both")
samples <- str_trim(unlist(str_split(opt$sample_id, pattern = ",")), 
                    side="both")
names(indir) <- samples
df <- tibble()
keep <- c()
for (n in names(indir)){
  tmp <- as_tibble(read.table(gzfile(file.path(indir[n], "ambient_genes.txt.gz")), 
                              header = TRUE, stringsAsFactors = FALSE))
  tmp$sample <- n
  df <- rbind(df, tmp)
  keep <- unique(c(keep,tmp %>% filter(top==TRUE) %>% 
                     pull(ensembl) %>% as.character()))
}
df <- df %>% filter(ensembl %in% keep)

#  Build matrix
flog.info("Creating ambient profile matrix")
df.wide <- df %>% dplyr::select(count_percentage, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = count_percentage) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide

# Save as ambient RNA profile
flog.info("Saving ambient rna profiles")
m.out = data.frame(gene=rownames(m), m)
write.table(m.out, file = file.path(opt$outdir, "ambient_rna_profile.tsv"),
            sep="\t", quote = FALSE, row.names = FALSE)

flog.info("Completed")
