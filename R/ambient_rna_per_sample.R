#' ---
#' title: "Ambient RNA"
#' output: 
#'  html_document:
#'   self_contained: false
#'   toc: true
#'   toc_float: true
#' params:
#'  task_yml: "ambient_rna.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "ambient_rna.log"
#' ---
#' ---
#' Run ambient RNA
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(yaml),
          require(DropletUtils),
          require(ggplot2),
          require(futile.logger),
          require(knitr),
          require(dplyr),
          require(ggrepel),
          require(ggExtra),
          require(cowplot),
          require(scales),
          require(reshape2)
          )

# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following parameters:
default_options <- list(
  # The path to the directory containing the cellranger 
  # raw matrix and associates barcodes and features
  "cellranger_dir" = "",
  
  # Path to the ambient RNA output file
  "outdir" = "ambient_rna.dir/",
  
  # Sample name
  "sample_name" = NULL,
  
  # Maximum UMI threshold to identify ambient droplets
  "umi" = 100,
  
  # Path to the file containing barcodes for blacklisting
  "blacklist" = NULL
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

flog.info("Reading cellranger input")
sce <- read10xCounts(file.path(opt$cellranger_dir))
flog.info("Dimensions of input object: %s", dim(sce))

## ######################################################################### ##
## ############## (ii) Filter data including blacklisting ################## ##
## ######################################################################### ##

if (!is.null(opt$blacklist)) {
  stopifnot(file.exists(opt$blacklist))
  flog.info("Read blacklisted cells from: %s", opt$blacklist)
  cells_blacklist <- data.frame(barcodes_in = scan(opt$blacklist, "character"))
  cells_blacklist$barcodes_blacklist = paste(cells_blacklist$barcodes_in, "1", sep="-")
  flog.info("Number of cells in given blacklist: %s", nrow(cells_blacklist))
  cells_remove <- colData(sce)$Barcode %in% cells_blacklist$barcodes_blacklist
  sce <- sce[, !cells_remove]
  flog.info("Number of cells after blacklisting: %s", ncol(sce))
} else {
  flog.info("No cells were blacklisted.") }

sce <- sce[,colSums(counts(sce))>0]
sce <- sce[rowSums(counts(sce))>0,]
flog.info("Dimensions of filtered object: %s", dim(sce))

my.counts <- counts(sce)

## ######################################################################### ##
## ################## (iii) Identify ambient droplets ###################### ##
## ######################################################################### ##

# Get cells with a maximum UMI count threshold
ambient.barcodes <- my.counts
ambient.barcodes <- DropletUtils:::.rounded_to_integer(ambient.barcodes)
max_umi = opt$umi
umi.sum <- as.integer(round(colSums(ambient.barcodes)))
ambient <- umi.sum <= max_umi
ambient.barcodes <- ambient.barcodes[, ambient]
flog.info("Number of droplets with less than %s umis: %s", 
          opt$umi, format(ncol(ambient.barcodes), big.mark = ","))
flog.info("This represents %s percent of the data", 
          round((ncol(ambient.barcodes)/ncol(my.counts)*100), digits=2))

n <- format(ncol(ambient.barcodes), big.mark = ",")
p <- round((ncol(ambient.barcodes)/ncol(my.counts)*100), digits=2)

#' ## Ambient droplets `r opt$sample_name`
#' The number of ambient droplets with less than `r max_umi` UMI counts is `r n`, which represents `r p` percent of the data.
#' 

## ######################################################################### ##
## ############### (iv) UMI distribution of ambient droplets ############### ##
## ######################################################################### ##

flog.info("Creating UMI count distribution plot for ambient droplets")
numis_per_cell <- colSums(ambient.barcodes)
numis_per_cell <- data.frame(nUMIs = numis_per_cell)

x <- ggplot(numis_per_cell, aes(x=nUMIs)) +
  geom_histogram(fill="white", col="black", bins=100) +
  scale_y_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::comma) +
  ylab("Barcode count") +
  xlab("UMI count") +
  theme_bw() +
  ggtitle("UMI count distribution of ambient RNA") + 
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        plot.title = element_text(size=17, hjust = 0.5, face="bold")) 

## ###################################################################### ##
## ############## (v) Genes detected in ambient droplets ############### ##
## ###################################################################### ##

flog.info("Identifying genes in ambient droplets")
# Get a gene dataframe with gene counts across ambient barcodes
ambient.genes <- tibble(ensembl=names(rowSums(ambient.barcodes)),
                            count=rowSums(ambient.barcodes))
# Get a gene dataframe with the number of ambient barcodes detecting a gene
nbarcodes_detected <- tibble(ensembl=names(rowSums(ambient.barcodes>0)),
                             n_barcodes_expressing=rowSums(ambient.barcodes>0))
# Join and get percentage of barcodes expressing a gene
ambient.genes <- left_join(ambient.genes, nbarcodes_detected, by="ensembl")
ambient.genes$percentage_barcodes_expressing <- (ambient.genes$n_barcodes_expressing/
                                                   ncol(ambient.barcodes))*100
# Keep only genes detected in ambient barcodes
ambient.genes <- ambient.genes[ambient.genes$count>0,]
# Get the percentage of counts attributed to each gene across all the ambient barcodes
ambient.genes$count_percentage <- ambient.genes$count/sum(ambient.genes$count)*100
# Add gene symbol information
ambient.genes <- left_join(ambient.genes, 
                           as.data.frame(rowData(sce)), 
                           by=c("ensembl"="ID"))
ambient.genes <- ambient.genes %>% dplyr::arrange(desc(count))
# Label the top 50 (or top 5% percent)
ambient.genes$label <- ""
top <- min(50,ceiling(nrow(ambient.genes)*5/100))
ambient.genes$label[1:top] <- ambient.genes$Symbol[1:top]
ambient.genes$top <- FALSE
ambient.genes$top[1:top] <- TRUE

flog.info("Plotting gene counts vs proportion of ambient droplets detecting a gene")
p <- ggplot(ambient.genes, aes(x=count, y=percentage_barcodes_expressing, label=label)) +
  geom_point(alpha=0.25) +
  geom_text_repel(segment.alpha =0.25) +
  scale_x_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::comma) +
  xlab("count") + 
  ylab("% of barcodes with 1 or more counts") + 
  theme_bw() +
  ggtitle("Gene UMI count in ambient droplets") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=16),
        plot.title = element_text(size=17, hjust = 0.5, face="bold")) 
p <- ggMarginal(p, type="density", margins = "x")

ambient.genes$Type <- NULL
ambient.genes$label <- NA
write.table(ambient.genes, gzfile(paste0(opt$outdir, "/ambient_genes.txt.gz")), 
            row.names = FALSE, quote = FALSE, sep="\t")

x <- cowplot::plot_grid(x,p, ncol=2, scale = 0.9, rel_widths = c(1,1))

#' ## UMIs in ambient droplets `r opt$sample_name`
#+ umi_in_ambient_droplets, include=TRUE, fig.width=12, fig.cap="", fig.align="center"
x
#'

# Barplot top 50 genes
flog.info("Creating barplots of top ambient RNA genes")
top_genes <- ambient.genes %>% dplyr::arrange(desc(count)) %>% head(50)
top_genes$type <- "Other"
top_genes$type[grep(paste0(c("^Rpl","^Rps"), collapse="|"), top_genes$Symbol, ignore.case=TRUE)] <- "Ribosomal"
top_genes$type[grep("^mt-", top_genes$Symbol, ignore.case = TRUE)] <- "Mitochondrial"
top_genes$Symbol <- factor(top_genes$Symbol, levels=top_genes$Symbol)


p1 <- ggplot(top_genes, aes(x=Symbol, y=count, fill = type)) +
  geom_bar( col="black",stat="identity") +
  scale_fill_manual(values = c("Other" = "grey", 
                               "Ribosomal" = "cadetblue2", 
                               "Mitochondrial" = "darkorange1")) + 
  theme_bw() +
  ylab("Total count") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=17),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = scales::comma)

p2 <- ggplot(top_genes, aes(x=Symbol, y=count_percentage, fill = type)) +
  geom_bar(col="black",stat="identity") +
  scale_fill_manual(values = c("Other" = "grey", 
                               "Ribosomal" = "cadetblue2", 
                               "Mitochondrial" = "darkorange1")) + 
  ylab("% of total \n UMI count") + 
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=17),
        legend.title = element_blank(),
        legend.position = "none") 

p3 <- ggplot(top_genes, aes(x=Symbol, y=percentage_barcodes_expressing, fill=type)) +
  geom_bar(col="black",stat="identity") +
  scale_fill_manual(values = c("Other" = "grey", 
                               "Ribosomal" = "cadetblue2", 
                               "Mitochondrial" = "darkorange1")) +
  theme_bw() +
  ylab("% of barcodes with \n 1 or more counts") + 
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=17),
        legend.title = element_blank(),
        legend.position = c(0.95, 0.75))

x <- cowplot::plot_grid(p1,p2,p3, nrow=3, align = "v")

#' ## Top genes detected in ambient RNA `r opt$sample_name`
#' The plots below show the top 50 genes with highest UMI counts across ambient droplets
#+ top_ambient_genes_barplots, include=TRUE, fig.width=15, fig.height=10, fig.cap="", fig.align="center"
x
#'

# Expression plot
flog.info("Creating boxplot of expression of top ambient RNA genes")
ambient.barcodes.topgenes <- ambient.barcodes[top_genes$ensembl,]
ambient.barcodes.topgenes <- reshape2::melt(as.matrix(ambient.barcodes.topgenes))
colnames(ambient.barcodes.topgenes) <- c("ensembl", "barcode", "value")
ambient.barcodes.topgenes <- left_join(ambient.barcodes.topgenes, 
                                       as.data.frame(rowData(sce)), 
                                       by=c("ensembl"="ID"))
ambient.barcodes.topgenes$type <- "Other"
ambient.barcodes.topgenes$type[grep(paste0(c("^Rpl","^Rps"), collapse="|"),
                                    ambient.barcodes.topgenes$Symbol, ignore.case = TRUE)] <- "Ribosomal"
ambient.barcodes.topgenes$type[grep("^mt-", 
                                    ambient.barcodes.topgenes$Symbol, ignore.case=TRUE)] <- "Mitochondrial"
ambient.barcodes.topgenes$Symbol <- factor(ambient.barcodes.topgenes$Symbol, 
                                           levels=unique(ambient.barcodes.topgenes$Symbol))

p <- ggplot(ambient.barcodes.topgenes, aes(x=Symbol, y=value, col=type)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::comma) +
  scale_color_manual(values = c("Other" = "grey", 
                               "Ribosomal" = "cadetblue2", 
                               "Mitochondrial" = "darkorange1")) +
  theme_bw() +
  ylab("Counts per\nbarcode") + 
  guides(colour = guide_legend(nrow = 1)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=17),
        legend.title = element_blank(),
        legend.position = c(0.87, 0.85))

#' ## Expression of top genes detected in ambient RNA `r opt$sample_name`
#+ top_ambient_genes_expression_boxplots, include=TRUE, fig.width=15, fig.height=4, fig.cap="", fig.align="center"
p
#'

flog.info("Completed")


