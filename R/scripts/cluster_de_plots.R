## Characterise a set of differentially expressed genes (per cluster)
## This script will produce:
## (1) MA plot (colored by p-value)
## (2) MA style plot for frequency of expression (colored by p-value)
## (3) Fold change vs Frequency change plot (coloured by p-value)
## (4) Traditional volcano plot
##
## The idea is to allow evaluation of the performance of the
## DE alogorithm applied (to aid/inform alogorithm selection)
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(optparse),
  #require(gplots),
  #require(reshape2),
  require(xtable),
  require(ggplot2),
  require(ggrepel),
  require(gridExtra),
  require(dplyr),
  require(cellhub)
)

# Options ----

option_list <- list(
    make_option(c("--degenes"), default="begin.Robj",
                help="Summary table of differentially expressed genes"),
    make_option(c("--cluster"), default="none",
                help="The cluster to characterise"),
    make_option(c("--testfactor"), default="NULL",
                help="a metadata factor used to group the violin plots"),
    make_option(c("--a"), default="NULL",
                help="first level of contrast"),
    make_option(c("--b"), default="NULL",
                help="second level of contrast"),
    make_option(c("--plotdirvar"), default="clusterMarkerDEPlotsDir",
                help="latex var containing location of the plots"),
    make_option(c("--pdf"), default=FALSE,
                help="Produce pdf plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$testfactor=="NULL") { opt$testfactor <- NULL }

cat("Running with options:\n")
print(opt)

## set up the ggplot theme
theme_set(theme_classic(base_size = 10))

if(is.null(opt$testfactor))
{
     file_suffix = NULL
#     a = "mean"
#     b = "other_mean"
 } else {
     file_suffix = paste0("between")
#     a = paste(opt$a,"mean",sep="_")
#     b = paste(opt$b,"mean",sep="_")
}

## initialise the text snippet
tex = ""

## read in the de genes
cluster <- opt$cluster
data <- read.table(gzfile(opt$degenes),header=T,as.is=T,sep="\t")

if(!cluster %in% data$cluster)
{
    ## we have no DE genes.
    message("No significant genes for cluster: ",cluster)
} else

{
data <- data[data$cluster == cluster,]

message("processing cluster: ",cluster)



## analyse expression level vs fold change


# note that the FC used here includes the specified pseudocount
# (see cluster_markers.py)

data$A <- 1/2 * (log2(data$mean_exprs) + log2(data$mean_exprs_other))

message("making expression MA")
print(head(data))
print(is.numeric(data$log2FC))
ma <- plotMA(data,xlab="expression level (log2)",
            m_col="log2FC", ylab="fold change (log2)")

message("making expression volcano")
vol <- plotVolcano(data,m_col="log2FC",
                   xlab="fold change (log2)",
                   ylab="adjusted p-value (-log10)")

## analyse percentage vs change in percentage

data$pctM <- log2((data$pct*100 + 1) / (data$pct_other*100 + 1))
data$A <- 1/2 * (log2(data$pct*100 + 1) + log2(data$pct_other*100 + 1))


message("making pct MA")
pma <- plotMA(data, m_col="pctM",
              xlab="percent cells (log2)",  ylab="percent change (log2)")

message("making pct volcano")
pvol <- plotVolcano(data, m_col="pctM", 
                    xlab="percent change (log2)",ylab="adjusted p-value (-log10)")

## directly compare fold change and percentage change
fve <- plotFvE(data, m_col="log2FC",freq_col="pctM",)

## sneakily extract the legend...
legend <- g_legend(fve)
fve <- fve + theme(legend.position = 'none')

gs <- list(ma, vol, pma, pvol, fve, legend)
rm(ma, vol, pma, pvol, fve, legend)

lay <- rbind(c(1,2),c(3,4),c(5,6))
ga <- arrangeGrob(grobs = gs, layout_matrix = lay)
rm(gs)

defn <- paste(c("dePlots",cluster,file_suffix),collapse=".")
defpath <- file.path(opt$outdir, defn)

save_ggplots(defpath,
             ga,
             to_pdf=opt$pdf,
             width=8,
             height=10)

rm(ga)
gc()

## start building figure latex...
subsectionTitle <- getSubsectionTex(paste0("Cluster ",cluster,": summary plots"))
tex <- c(tex, subsectionTitle)

deCaption <- paste("Differential expression summary plots for cluster ",cluster)
tex <- c(tex, getFigureTex(defn, deCaption,plot_dir_var=opt$plotdirvar))

# tex <- c(tex, pos_tex, neg_tex)

tex_file <- file.path(opt$outdir,
                      paste(c("characterise.degenes",cluster,file_suffix,"tex"),collapse="."))

writeTex(tex_file, tex)

}

message("completed")
