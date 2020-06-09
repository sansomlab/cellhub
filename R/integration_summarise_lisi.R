#' ---
#' title: "Summary of iLISI calculations across tools"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: ""
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
stopifnot(require(knitr))
stopifnot(require(RColorBrewer))
stopifnot(require(futile.logger))

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

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))

flog.info("Running with options:", opt, capture = TRUE)
flog.info("\n")


## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

plot_data <- data.frame()
for (f in strsplit(opt$lisi_files, ",")[[1]]){
  df = read.csv(f)
  df = df[, !colnames(df) %in% c("UMAP_1", "UMAP_2", "LARGEVIS_1", "LARGEVIS_2")]
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
                                        levels=c(unique(plot_data$integration_setting[grepl("rawdata",plot_data$integration_setting)]), 
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
# here the first 8 colors of the palette are selected and then potentially spread
# if more levels exist
colors_use <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(length(levels))



gp <- ggplot(plot_data, aes(x=lisi_value, color = integration_setting)) + geom_density() 
gp <- gp + theme_bw() + scale_color_manual(values = colors_use) 
gp <- gp + xlab("iLISI")

ggsave(file.path(opt$outdir, "iLISI_all.pdf"), 
       plot = gp, device = cairo_pdf, width=15,
       units="in")


#' 
#' **Summary of iLISI results:**
#+ summary_lisi, include=TRUE, fig.height=3, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE



gp <- ggplot(plot_data, aes(x=lisi_value, color = parameters)) + geom_density()
gp <- gp + facet_wrap(~tool, ncol=1) + scale_color_manual(values = colors_use) + theme_bw()
gp <- gp + xlab("iLISI")

ggsave(file.path(opt$outdir, "iLISI_by_tool.pdf"), 
       plot = gp, device = cairo_pdf, width=10,
       height=10, units="in")


#'  
#'  **Summary of iLISI results split by integration tool:**
#+ summary_lisi_density_byTool, include=TRUE, fig.height=5, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE


gp <- ggplot(plot_data, aes(x=lisi_value, fill = parameters)) + geom_histogram(position = "dodge")
gp <- gp + facet_wrap(~tool, ncol=1) + scale_fill_manual(values = colors_use) + theme_bw()
gp <- gp + xlab("iLISI") + ylab("Number of cells")

ggsave(file.path(opt$outdir, "iLISI_histograms_by_tool.pdf"), 
       plot = gp, device = cairo_pdf, width=10,
       height=10, units="in")

#' 
#'  **Summary of iLISI results split by integration tool (histogram):**
#+ summary_lisi_histogram_byTool, include=TRUE, fig.height=5, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE


write.table(plot_data, file.path(opt$outdir, "iLISI_all.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

flog.info("Completed")


