#' ---
#' title: "Run kBET"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "kbet.log"
#' ---
#' ---
#' Run kBET post-integration
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(Seurat))
stopifnot(require(Matrix))
stopifnot(require(dplyr))
stopifnot(require(tenxutils))
stopifnot(require(knitr))
stopifnot(require(futile.logger))
stopifnot(require(kBET))

# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Name of folders to output rds, pdf and png files
  "outdir" = "",
  
  # Path to Seurat object
  "seurat_obj" = "",
  
  # Variable for splitting seurat object for integration
  "split_var" = "sample_id",
  
  # Correct assay to use
  "assay" = "integrated",
  
  # Correct dimension reduction slot to use
  "slot" = "pca",
  
  # Number of PCs to use, kbet uses 50 by default, here 30 used
  "nPC" = 30
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

s_integrated <- readRDS(file = opt$seurat_obj)

# write out info about assay and slot that were selected
flog.info("Chosen assay is %s", opt$assay)
flog.info("Chosen slot is %s", opt$slot)

# change the default assay if required as MixingMetric function does not have
# an assay argument
if (DefaultAssay(s_integrated) != opt$assay){
  flog.info("Reset the DefaultAssay because it did not match the 
            required assay.")
  DefaultAssay(s_integrated) <- opt$assay
}

## ######################################################################### ##
## ###################### (ii) run kBET #################################### ##
## ######################################################################### ##

#' ## Assessment of integration using kBET
#' 
#' Paper for kBET can be found [here](https://www.nature.com/articles/s41592-018-0254-1)

flog.info("Running kbet using variable: %s", opt$split_var)

data_kbet = Embeddings(s_integrated, reduction = opt$slot)

if (ncol(data_kbet) < opt$nPC) {
    opt$nPC <- ncol(data_kbet)
    flog.info("Change number of components as not sufficient number present ...")
}

data_kbet = data_kbet[,1:opt$nPC]
metadata = s_integrated[[]][,opt$split_var]

batch.estimate <- kBET(data_kbet, metadata, plot=FALSE, do.pca = FALSE)

## ######################################################################### ##
## ###################### (ii) make plots ################################## ##
## ######################################################################### ##


flog.info("Make kbet plots")

## make the standard kbet plot
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))



gp <- ggplot(plot.data, aes(class, data)) + geom_boxplot() 
gp <- gp + labs(x='Test', y='Rejection rate',title='kBET test results') 
gp <- gp + theme_bw() + scale_y_continuous(limits=c(0,1))

save_ggplots(file.path(opt$outdir, "k_bet_boxplot"), 
             gp)

#' Boxplot of kBET results
#+ boxplot_kbet, include=TRUE, fig.height=3
print(gp)
#+ include=FALSE
#'


flog.info("Completed")
