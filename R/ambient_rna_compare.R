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

stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(futile.logger))
stopifnot(require(knitr))
stopifnot(require(dplyr))
stopifnot(require(ggrepel))
stopifnot(require(stringr))
stopifnot(require(tidyverse))
stopifnot(require(ComplexHeatmap))
stopifnot(require(RColorBrewer))
stopifnot(require(grDevices))
stopifnot(require(optparse))


# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path)

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
  "sample_table" = "",
  
  # Additional columns from the sample table (input_samples.txt) to show in summary heatmaps
  # Requires comma-separated column names from the input_samples.txt to specify the annotations to plot
  "plot_annotation" = ""
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

## ################################################################ ##
## ###################### (ii) Make heatmaps ###################### ##
## ################################################################ ##

# Prepare annotation for heatmap
sample_table <- read.table(opt$sample_table, header=TRUE, stringsAsFactors = FALSE)
annotation <- sample_table %>% dplyr::select(-path)
annotation <- apply(annotation, 2, as.character) %>% as.data.frame()
rownames(annotation) <- annotation$sample_id
annotation$sample_id <- NULL
# Keep annotations from options
keep <-  str_trim(unlist(str_split(opt$plot_annotation, 
                                     pattern = "," )),
                  side="both")
annotation <- annotation[ , keep , drop=FALSE]
# Keep annotations with more than 1 level
annotation <- annotation[,which(apply(annotation,2, function(x) length(unique(x)))>1), drop=FALSE]

#  Make heatmap - percentage of total UMI count ---------------
flog.info("Creating heatmap of top ambient gene UMI count")
df.wide <- df %>% dplyr::select(count_percentage, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = count_percentage) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide

# Set NAs to 0
m.scaled <- m
m.scaled[is.na(m.scaled)] <- 0

# Row scale matrix
m.scaled <- t(scale(t(m.scaled)))

# Cluster genes 
d <- dist(m.scaled)
h <- hclust(d)
ro <- rownames(m.scaled)[h$order]

# Cluster samples
d <- dist(t(m.scaled))
h <- hclust(d)
co <- colnames(m.scaled)[h$order]

# Order matrix based on clustering
m.scaled <- m.scaled[ro,co]

# Make annotation in sample col order as matrix
if ( ncol(annotation)>0) {
  annotation <- annotation[co, ,drop = FALSE]
  # Color palettes
  palettes <- c("Paired", "Set2", "Dark2", "Set1", 
                "Accent", "Set3", "Pastel1", "Pastel2")
  col.list <- list()
  if (ncol(annotation) > 1) {
    for (c in 1:ncol(annotation)){
      levels <- unique(as.character(annotation[,c]))
      # here the first 8 colors of the palette are selected and then potentially spread
      # if more levels exist
      cols <- colorRampPalette(RColorBrewer::brewer.pal(palettes[c], n = 8))(length(levels))
      names(cols) <- levels
      col.list[[colnames(annotation)[c]]] <- cols
      }
    } else {
      levels <- unique(as.character(annotation[,1]))
      cols <- colorRampPalette(RColorBrewer::brewer.pal(palettes[1], n = 8))(length(levels))
      names(cols) <- levels
      col.list[[colnames(annotation)]] <- cols
    }
  ha <- HeatmapAnnotation(df = annotation, 
                          col = col.list)
} else {
  ha <- NULL
}

# Calculate size of plot
ws <- setNames(seq(from = 2.5, to=10, length.out = length(15:2)), nm = as.character(15:2))
w <- ws[as.character(length(samples))]

# Make heatmap
m1 <- Heatmap(m.scaled,
        name = "% of total UMI count\n(row-scaled)", 
        top_annotation = ha,
        column_names_rot = 65,
        row_names_gp = gpar(fontsize = 10),
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        column_title = "Percentage of\ntotal UMI count\n(row-scaled)",
        column_names_side = "top",
        rect_gp = gpar(col='white', lwd=0.5), 
        heatmap_width = unit(w, "npc"),
        na_col = "black")

#  Make heatmap - percentage of total UMI count, non-scaled ---------------
flog.info("Creating heatmap of top ambient gene barcode percentage detection")
df.wide <- df %>% dplyr::select(count_percentage, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = count_percentage) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide

# Order matrix based on clustering (percentage of total UMI count heatmap)
m <- m[ro,co]

# Plot
ws <- setNames(seq(from = 2.5, to=10, length.out = length(15:2)), nm = as.character(15:2))
w <- ws[as.character(length(samples))]

col_fun = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd")

m1ns <- Heatmap(m,
              name = "% of total UMI count ", 
              top_annotation = ha,
              column_names_rot = 65,
              row_names_gp = gpar(fontsize = 10),
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              column_title = "Percentage of\ntotal UMI count",
              column_names_side = "top",
              rect_gp = gpar(col='white', lwd=0.5),
              heatmap_width = unit(w, "npc"),
              col = col_fun, 
              na_col = "black")

# Save as ambient RNA profile
m.out = data.frame(gene=rownames(m), m)
write.table(m.out, file = file.path(opt$outdir, "ambient_rna_profile.tsv"),
            sep="\t", quote = FALSE, row.names = FALSE)

# Make heatmap - percentage_barcodes_expressing ---------------
flog.info("Creating heatmap of top ambient gene barcode percentage detection")
df.wide <- df %>% dplyr::select(percentage_barcodes_expressing, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = percentage_barcodes_expressing) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide

# Set NAs to 0
m.scaled <- m
m.scaled[is.na(m.scaled)] <- 0

# Row scale matrix
m.scaled <- t(scale(t(m.scaled)))

# Order matrix based on clustering (percentage of total UMI count heatmap)
m.scaled <- m.scaled[ro,co]

ws <- setNames(seq(from = 2.5, to=10, length.out = length(15:2)), nm = as.character(15:2))
w <- ws[as.character(length(samples))]

m2 <- Heatmap(m.scaled,
             name = "% of barcodes\nwith 1 or more counts\n(row-scaled)", 
             top_annotation = ha,
             column_names_rot = 65,
             row_names_gp = gpar(fontsize = 10),
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             column_title = "Percentage of\nbarcodes with 1 or more\ncounts (row-scaled)",
             column_names_side = "top",
             rect_gp = gpar(col='white', lwd=0.5),
             heatmap_width = unit(w, "npc"),
             na_col = "black")

# Make heatmap - percentage_barcodes_expressing, non-scaled ---------------
flog.info("Creating heatmap of top ambient gene barcode percentage detection")
df.wide <- df %>% dplyr::select(percentage_barcodes_expressing, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = percentage_barcodes_expressing) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide
m <- m[ro,co]

ws <- setNames(seq(from = 2.5, to=10, length.out = length(15:2)), nm = as.character(15:2))
w <- ws[as.character(length(samples))]

col_fun = RColorBrewer::brewer.pal(n = 8, name = "YlOrRd")

m2ns <- Heatmap(m,
              name = "% of barcodes\nwith 1 or more counts", 
              top_annotation = ha,
              column_names_rot = 65,
              row_names_gp = gpar(fontsize = 10),
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              column_title = "Percentage of barcodes\nwith 1 or more counts",
              column_names_side = "top",
              rect_gp = gpar(col='white', lwd=0.5),
              heatmap_width = unit(w, "npc"),
              col = col_fun,
              na_col = "black")


# Make heatmap - top or not ---------------
flog.info("Creating heatmap of top ambient genes in samples")
df.wide <- df %>% dplyr::select(top, Symbol, sample) %>% 
  pivot_wider(.,  names_from = sample,values_from = top) 
rownames.df.wide <- as.character(df.wide$Symbol)
df.wide$Symbol <- NULL
m <- as.matrix(df.wide)
rownames(m) <- rownames.df.wide
m[m==TRUE] = "True"
m[m==FALSE] = "False"
m <- m[ro,co]

ws <- setNames(seq(from = 2.5, to=10, length.out = length(15:2)), nm = as.character(15:2))
w <- ws[as.character(length(samples))]

m3 <- Heatmap(m,
              name = "Top ambient gene", 
              top_annotation = ha,
              column_names_rot = 65,
              row_names_gp = gpar(fontsize = 10),
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              column_title = "Top ambient gene",
              col = c("True" = "grey30", "False" = "grey90"),
              column_names_side = "top",
              rect_gp = gpar(col='white', lwd=0.5),
              heatmap_width = unit(w, "npc"),
              na_col = "black")


# Join the three heatmaps together, Percentage of total UMI count of top ambient genes ---------------
flog.info("Joinning heatmaps")
ht_list <- m1 + m1ns + m3
w = 110 * length(samples)

png(paste0(opt$outdir, "/top_ambient_heatmap.png"), height = 1000, width=w, type="cairo")
draw(ht_list, ht_gap = unit(0.5, "cm"))
dev.off()

# set chunk options
w = 14
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path,
               fig.width = w)


#' ## Percentage of total UMI count of top ambient genes
#' The percentage of UMI counts assigned to each gene in the ambient droplets was calculated in each sample. 
#' The top 50 ambient RNA genes (i. e. the genes with highest percentage of UMI count) were identifyied in each sample.
#' The heatmaps below show the set of genes that correspond to the union of the top 50 ambient RNA genes across samples. 
#' The heatmap on the left shows the percentage of UMI counts for the genes (value is scaled across samples).
#' The heatmap in the center shows the percentage of UMI counts for the genes (value is real percentage).
#' The heatmap on the right shows whether a gene was in the top 50 or not in each sample (values are True or False).
#' NA values are in black, and they mean that the gene was not present in the ambient droplets of the corresponding sample.
#+ top_ambient_genes_heatmap_p_total_umi, include=TRUE,  fig.height=15, fig.cap="", fig.align="center"
draw(ht_list, ht_gap = unit(50, "mm"))
#'

# Join the three heatmaps together, Percentage of barcodes with 1 or more counts of top ambient genes ---------------
ht_list <- m2 + m2ns + m3

#' ## Percentage of barcodes with 1 or more counts of top ambient genes
#' The percentage of barcodes expressing a gene (> 1 UMI count) was calculated in each sample.
#' The heatmaps below show the set of genes that correspond to the union of the top 50 ambient RNA genes across samples. 
#' The heatmap on the left shows the percentage of barcodes expressing each gene on the rows (value is scaled across samples).
#' The heatmap in the center shows the percentage of barcodes expressing each gene on the rows (value is real percentage).
#' The heatmap on the right shows whether a gene was in the top 50 or not in each sample (values are TRUE or FALSE).
#' NA values are in black, and they mean that the gene wasn't present in the ambient droplets of the corresponding sample.
#+ top_ambient_genes_heatmap_p_exprs_bc, include=TRUE,  fig.height=15, fig.cap="", fig.align="center"
draw(ht_list, ht_gap = unit(50, "mm"))
#'

# Other plots
# dirs <- paste0(indir)
# current_dirc <- getwd()
# for (d in dirs){
#   system(paste0("cd ", opt$outdir, "/fig.dir"))
#   system(paste0("ln -s ", current_dirc, 
#                 "/", d, "/ambient_umi_distribution.png ", 
#                 paste0(current_dirc,"/",opt$outdir), 
#                 "/fig.dir/", 
#                 gsub("\\.dir.*", "", d), 
#                 "_ambient_umi_distribution.png"))
#   system(paste0("cd ", current_dirc))
# }
# 
# inplots=paste0("fig.dir/", names(indir),"_ambient_umi_distribution.png")
# names(inplots) <- names(indir)
# for (p in inplots){
#   cat("![](",p,")")
#   cat("\n")
# }


flog.info("Completed")


