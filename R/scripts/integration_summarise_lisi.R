#' ---
#' title: "Summary of iLISI calculations across tools"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "summarise_lisi.log"
#' ---
#' ---
#' Summarise results from iLISI
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(tidyr))
stopifnot(require(futile.logger))
stopifnot(require(colormap))

# set chunk options
# opts_chunk$set(echo=FALSE,
#                warning=FALSE,
#                message = FALSE,
#                include = FALSE,
#                fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------

# these parameters are passed from Rscript run
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

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))
flog.info("Running with parameters: ", params, capture = TRUE)


# The script expects the following paramters from yml:
default_options <- list(
  # Name of folders to output rds, pdf and png files
  "outdir" = "",
  
  # Path to lisi output files
  "lisi_files" = ""
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

flog.info("Running with options:", opt, capture = TRUE)
flog.info("\n")


## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

plot_data <- data.frame()
for (f in strsplit(opt$lisi_files, ",")[[1]]){
  df = read.csv(f, sep="\t")
  df$file = f
  df = df %>% separate(file, c("sample", "integration_tool", "parameters", "assess_folder", "file_name"), sep = "/") %>%
        select (-sample, -assess_folder, -file_name) %>%
        separate(integration_tool, c("tool", "int", "dir"), sep = "\\.") %>% select(-int, -dir) %>%
        separate(parameters, c("parameters", "dir"), sep = "\\.run") %>% select(-dir) %>%
        unite("integration_setting", c(tool, parameters), remove = FALSE)
  plot_data = rbind(plot_data, df)
}

# reorder data for plotting
plot_data$tool <- factor(plot_data$tool, levels=c("rawdata", 
                                                      unique(plot_data$tool[!plot_data$tool == "rawdata"])))
plot_data$integration_setting <- factor(plot_data$integration_setting, 
                                        levels=c(unique(plot_data$integration_setting[grepl("rawdata", 
                                                                                             plot_data$integration_setting)]), 
                                                 unique(plot_data$integration_setting[!grepl("rawdata", 
                                                                                             plot_data$integration_setting)])))

## ######################################################################### ##
## ###################### (ii) Make plots #############@#################### ##
## ######################################################################### ##

#' ## Summary of iLISI results across different integration tools
#' 
#' The maximum iLISI is the number of levels the variables used for integration has, 
#' e.g. if there were two batches, the maximum iLISI is 2.     
#' 

levels <- unique(as.character(plot_data$integration_setting))
cm_samples <- colormap(colormap = colormaps$portland,
                       nshade = length(levels))
cm_parameters <- colormap(colormap = colormaps$portland,
                       nshade = length(unique(plot_data$parameters)))


lisi_types = colnames(plot_data)[grepl("iLISI", colnames(plot_data))]
flog.info("Lisi was run on variable: %s", paste(lisi_types, collapse=","))

for (l in lisi_types){
  flog.info("Plotting lisi summary for %s", l)
  gp <- ggplot(plot_data, aes_string(x=l, color = "integration_setting")) + geom_density() 
  gp <- gp + theme_bw() + scale_color_manual(values = cm_samples) 
  gp <- gp + xlab("iLISI") + ggtitle(l)
  
  ggsave(file.path(opt$outdir, paste0(l, ".pdf")),
         plot = gp, device = cairo_pdf,width = 10)
}

for (l in lisi_types){
  gp <- ggplot(plot_data, aes_string(x=l, color = "parameters")) + geom_density()
  gp <- gp + facet_wrap(~tool, ncol=1) + scale_color_manual(values = cm_parameters) + theme_bw()
  gp <- gp + xlab("iLISI") + ggtitle(l)
  
  ggsave(file.path(opt$outdir, paste0(l, "_byTool.pdf")),
         plot = gp, device = cairo_pdf, width = 8)
  }



for (l in lisi_types){
  gp <- ggplot(plot_data, aes_string(x=l, fill = "parameters")) + geom_histogram(position = "dodge")
  gp <- gp + facet_wrap(~tool, ncol=1) + scale_fill_manual(values = cm_parameters) + theme_bw()
  gp <- gp + xlab("iLISI") + ylab("Number of cells") + ggtitle(l)
  
  ggsave(file.path(opt$outdir, paste0(l, "_histogram_byTool.pdf")),
         plot = gp, device = cairo_pdf, width = 8, height = 8)
}


write.table(plot_data, file.path(opt$outdir, "iLISI_all.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

flog.info("Completed")


