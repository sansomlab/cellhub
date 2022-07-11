## Characterise a set of differentially expressed genes (per cluster)
## This script will produce:
## 
## Violin plots of the top DE genes.
##
## The idea is to allow evaluation of the performance of the
## DE alogorithm applied (to aid/inform alogorithm selection)
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(optparse),
  require(gplots),
  require(reshape2),
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
    make_option(c("--loom"), default="data.loom",
                help="path to the loom file"),
    make_option(c("--matrix_loc"), default="matrix",
                help="location of the matrix in the loom file"),
    make_option(c("--scale"), default=FALSE,
                help="should the expression data be scaled for the heatmap"),
    make_option(c("--barcode_id_loc="), default="col_attrs/barcode_id",
                help="location of the barcode ids in the loom file"),
    make_option(c("--gene_id_loc="), default="row_attrs/gene_ids",
                help="location of the gene ids in the loom file"),
    make_option(c("--clusterids"), default="none",
                help="A tsv file containing the cluster identities"),
    make_option(c("--metadata"), default="none",
                help="A tsv file containing the metadata"), 
    make_option(c("--cluster"), default="none",
                help="The cluster to characterise"),
    make_option(c("--testfactor"), default="NULL",
                help="a metadata factor used to group the violin plots"),
    make_option(c("--a"), default="NULL",
                help="first level of constrast"),
    make_option(c("--b"), default="NULL",
                help="second level of constrast"),
    make_option(c("--plotdirvar"), default="clusterMarkerDEPlotsDir",
                help="latex var containing location of the plots"),
    make_option(c("--useminfc"), default=FALSE,
                help="Use minimum fold change for second page of violin plots"),
    make_option(c("--ncol"), type="integer", default=4,
                help="Number of columns for the plots"),
        make_option(c("--nrow"), type="integer", default=4,
                help="Number of rows for the plots"),
    make_option(c("--pointsize"), type="integer", default=FALSE,
                help="pointsize for the violin plots"),
    make_option(c("--pdf"), default=FALSE,
                help="Produce pdf plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$testfactor=="NULL") { opt$testfactor <- NULL }

cat("Running with options:\n")
print(opt)

# set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")

## set up the ggplot theme
theme_set(theme_classic(base_size = 10))

if(is.null(opt$testfactor))
{
    file_suffix = NULL
    a = "cluster_mean"
    b = "other_mean"
} else {
    file_suffix = paste0("between")
    a = paste(opt$a,"mean",sep="_")
    b = paste(opt$b,"mean",sep="_")
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


## start building figure latex...
subsectionTitle <- getSubsectionTex(paste0("Cluster ",cluster,": summary plots"))
tex <- c(tex, subsectionTitle)

deCaption <- paste("Differential expression summary plots for cluster ",cluster)
tex <- c(tex, getFigureTex(defn, deCaption,plot_dir_var=opt$plotdirvar))

## make two pages of violin plots for the top average cluster markers:
## Page (1) (a) top +ve by p-value
##          (b) top +ve by fold change (non redundant with (a))
## Page (2) (a) top -ve by p-value
##          (b) top -ve by fold change (non redundant with (b))
##
## - Only significant genes will be plotted
## - For cluster marker genes, min fold change (vs all other clusters is used).

## order by should be one of "p-value", "fold-change" or "min-fold-change".
## enforce significance

## read in the seurat object
cluster_ids <- read.table(opt$clusterids, header=T, sep="\t")

violin_fn_template <- paste(c("violinPlots.type",cluster,file_suffix),collapse=".")

data <- data[data$p.adj < 0.05,]

if(opt$useminfc)
{   fc_type="minimum fold-change vs all other clusters"
} else {    fc_type="fold change" }

ncol = opt$ncol

if(is.null(opt$testfactor))
{
    analysis_title = " marker genes"
    ident.include <- NULL
} else {
    analysis_title = "ly differentially expressed genes"
    ident.include <- cluster
}

## make the +ve plots
message("making violin plots for +ve genes")
pos_tex <- violinPlotSection(data=data, 
                             loom=opt$loom,
                             matrix_loc=opt$matrix_loc,
                             barcode_id_loc=opt$barcode_id_loc,
                             scale=opt$scale,
                             gene_id_loc=opt$gene_id_loc,
                             cluster_ids=opt$clusterids, 
                             type="positive",
                             group.by = opt$testfactor,
                             ident.include = ident.include, vncol=opt$ncol, vnrow=opt$nrow,
                             pt_size=opt$pointsize,
                             outdir = opt$outdir,
                             analysis_title = analysis_title, fc_type = fc_type,
                             use.minfc = opt$useminfc,
                             plot_dir_var=opt$plotdirvar,
                             to_pdf=opt$pdf)

## make the -ve plots
message("making violin plots for -ve genes")
neg_tex <- violinPlotSection(data=data, 
                             loom=opt$loom,
                             matrix_loc=opt$matrix_loc,
                             barcode_id_loc=opt$barcode_id_loc,
                             scale=opt$scale,
                             gene_id_loc=opt$gene_id_loc,
                             cluster_ids=opt$clusterids,
                             type="negative",
                             group.by = opt$testfactor,
                             ident.include = ident.include, vncol=opt$ncol, vnrow=opt$nrow,
                             pt_size=opt$pointsize,
                             outdir = opt$outdir,
                             analysis_title = analysis_title, fc_type = fc_type,
                             use.minfc = opt$useminfc,
                             plot_dir_var=opt$plotdirvar,
                             to_pdf=opt$pdf)

tex <- c(tex, pos_tex, neg_tex)

tex_file <- file.path(opt$outdir,
                      paste(c("characterise.degenes",cluster,file_suffix,"tex"),collapse="."))

writeTex(tex_file, tex)

}


message("completed")
