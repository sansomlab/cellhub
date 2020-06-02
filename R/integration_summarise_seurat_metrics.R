#' ---
#' title: "Summary of Seurat metrics across tools"
#' output:
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/kjansen/integration_pipeline/Rmd/integration_harmony.test.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "summarise_seurat_metrics.log"
#' ---
#' ---
#' Summarise Seurat metrics from the different integration tools
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
stopifnot(require(reshape2))
stopifnot(require(tenxutils))
stopifnot(require(RColorBrewer))
stopifnot(require(grid))
stopifnot(require(gridExtra))
stopifnot(require(knitr))
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

  # Path to metrics output files
  "metrics_files" = "",

  # Variable for splitting seurat object for integration
  "split_var" = "sample_id"
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

plot_data <- data.frame()
for (f in strsplit(opt$metrics_files, ",")[[1]]){
    flog.info("Adding plot data for metrics file:", f, capture=TRUE)
  df = read.csv(f, stringsAsFactors = FALSE)
  df$file = f
  df = df %>% separate(file, c("sample", "integration_tool", "parameters", "assess_folder", "file_name"), sep = "/") %>%
    select (-sample, -assess_folder, -file_name) %>%
    separate(integration_tool, c("tool", "int", "dir"), sep = "\\.") %>% select(-int, -dir) %>%
    separate(parameters, c("parameters", "dir"), sep = "\\.run") %>% select(-dir) %>%
    unite("integration_setting", c(tool, parameters), remove = FALSE)
  plot_data = rbind(plot_data, df)
}

flog.info("head of plot data:", head(plot_data), capture=TRUE)

# reorder data for plotting
plot_data$method <- factor(plot_data$method, levels=c("rawdata",
                                                      unique(plot_data$method[!plot_data$method == "rawdata"])))

plot_data$integration_setting <- factor(plot_data$integration_setting,
                                        levels=c(unique(plot_data$integration_setting[grepl("rawdata",
                                                                                             plot_data$integration_setting)]),
                                                 unique(plot_data$integration_setting[!grepl("rawdata",
                                                                                            plot_data$integration_setting)])))

## ######################################################################### ##
## ###################### (ii) Make plots ################################## ##
## ######################################################################### ##

flog.info("Make plots ...")

#' ## Summary of Seurat metrics from all tools
#' Please note that the local structure cannot be determined for the rawdata and
#' is therefore set to 0.\
#'
#' **Seurat Mixing Metric for each integration variable individually:**

levels <- unique(as.character(plot_data$integration_setting))
# here the first 8 colors of the palette are selected and then potentially spread
# if more levels exist
colors_use <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(length(levels))


# plot each integration grouping individually [only use shape if less than 10 samples]
if (length(unique(plot_data[,opt$split_var])) < 11){
  gp <- ggplot(data=plot_data, aes_string(x="mean_mixing_metric",
                                        y="meanLocalStructure", shape = opt$split_var,
                                        color = "integration_setting"))
  gp <- gp + scale_shape_manual(values = c(15:18, 1:14))
} else {
  gp <- ggplot(data=plot_data, aes_string(x="mean_mixing_metric",
                                          y="meanLocalStructure", 
                                          color = "integration_setting"))
}
gp <- gp + geom_point(size=1) + scale_color_manual(values = colors_use)
gp <- gp  + xlab("Mean Mixing Metric for each condition")
gp <- gp +theme_bw() 

save_ggplots(file.path(opt$outdir, "Alignment_metrics_summary_by_grouping_var"),
             gp, width = 10)



#'
#+ summary_seurat_byVariable, include=TRUE, fig.height=3, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE

write.csv(file.path(opt$outdir, "summary_alignment_metrics_seurat_package.csv"))

#' **Seurat Mixing Metric for each integration tool (mean across variables):**
#'

# mean across each tool
plot_data = plot_data %>% group_by(integration_setting,method,parameters) %>% summarise(mean_MM = mean(mean_mixing_metric),
                                                       mean_LS = mean(meanLocalStructure))


gp <- ggplot(data=plot_data, aes(x=mean_MM, y=mean_LS, color = integration_setting)) + geom_point(size=2)
gp <- gp + coord_cartesian(ylim = c(0,1)) + theme_bw() + scale_color_manual(values = colors_use)
gp <- gp + xlab("Mean Mixing Metric for each tool") + ylab("Mean Local Structure")

save_ggplots(file.path(opt$outdir, "Alignment_metrics_summary"),
             gp, width = 10)


#+ summary_seurat, include=TRUE, fig.height=4, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE

gp <- ggplot(data=plot_data, aes(x=mean_MM, y=mean_LS, color = parameters)) + geom_point(size=2)
gp <- gp + facet_wrap(~method) + theme_bw() + scale_color_manual(values = colors_use)
gp <- gp + xlab("Mean Mixing Metric for each tool") + ylab("Mean Local Structure")
gp <- gp + coord_cartesian(ylim = c(0,1))

save_ggplots(file.path(opt$outdir, "Alignment_metrics_summary_by_tool"),
             gp)

#' **Seurat Mixing Metric facetted by integration tool (mean across variables):**
#'
#+ summary_seurat_byTool, include=TRUE, fig.height=4, fig.width=8, fig.cap="", fig.align="center"
print(gp)
#+ include=FALSE
#'

write.table(plot_data, file.path(opt$outdir, "seurat_metrics_all.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

flog.info("Completed")
