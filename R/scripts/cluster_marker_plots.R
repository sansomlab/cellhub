stopifnot(
  require(optparse),
  require(cellhub),
  require(ggplot2),
  require(reshape2),
  require(dplyr),
  require(colormap),
  require(circlize),
  require(ComplexHeatmap)
)


# Options ----

option_list <- list(
    make_option(c("--markers"),default="none",
                help="A table containing the marker information"),
    make_option(c("--loom"), default="data.loom",
                help="path to the loom file"),
    make_option(c("--scaled_matrix_loc"), default="matrix",
                help="location of the scaled matrix in the loom file"),
    make_option(c("--data_matrix_loc"), default="layers/log1p",
                help="location of the normalised data matrix in the loom file"),
    make_option(c("--scale"), default=FALSE,
                help="should the expression data be scaled for the heatmap"),
    make_option(c("--barcode_id_loc="), default="col_attrs/barcode_id",
                help="location of the barcode ids in the loom file"),
    make_option(c("--gene_id_loc="), default="row_attrs/gene_name",
                help="location of the gene ids in the loom file"),
    make_option(c("--clusterids"), default="none",
                help="A tsv file containing the cluster identities"),
    make_option(c("--metadata"), default="none",
                help="A tsv file containing the metadata"),
    make_option(c("--annotation"), default="none",
                help="A tsv file containing the annation"),
    make_option(c("--group"), default=NULL,
                help="A column in the cell metadata used to define the groups"),
  make_option(c("--cluster"), default=0,
                help="The cluster to generate the plots for"),
  make_option(c("--rdimstable"), default=NULL,
                help="the reduced dimension coordinates"),
  make_option(c("--rdim1"), default="UMAP_1", 
                help="Rdims x coordinate column"),
  make_option(c("--rdim2"), default="UMAP_2",
                help='Rdims y coordinate'),
  make_option(c("--ngenes"), default=24, type="integer",
                help='number of genes to pull, currently ignored...'),
  make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir"),
  make_option(c("--pdf"), default=FALSE,
                help="outdir"))

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

outprefix = file.path(opt$outdir,paste("cluster",opt$cluster,sep="."))

message("reading in the markers")
x <- read.table(opt$markers, sep="\t", header=T, as.is=T)

message("reading in the rdims table")
rdims <- read.table(opt$rdimstable,
                   sep="\t", header=T,  as.is = T)

message("filtering the markers")
x <- x[x$cluster==as.numeric(opt$cluster) & x$p.adj<0.1 & !is.na(x$p.adj),]

#print(head(x))

n_select = 16

# pull out by:
# (i) top p-value vs all other (ranking by the test statistic
#     score to avoid ties)
# (ii) and minimum log fold change against each other cluster
top_by_pval <- x %>% arrange(desc(score)) %>%
                   slice_head(n = n_select)

n_select = 8
top_by_min_log2FC <- x[!x$gene_name %in% top_by_pval$gene_name,] %>%
                     arrange(desc(min_log2FC)) %>%
                     slice_head(n = n_select)

x <- x[x$gene_name %in% c(top_by_pval$gene_name,
                     top_by_min_log2FC$gene_name),]

# marker gene heatmap.

if(nrow(x) > 0) {

print(head(x,n=6))

message("making the heatmap")
mch <- markerComplexHeatmap(loom_path=opt$loom,
                            matrix_loc=opt$scaled_matrix_loc,
                            barcode_id_loc=opt$barcode_id_loc, #"col_attrs/barcode_id",
                            # gene_id_loc=opt$gene_id_loc, # "row_attrs/gene_name",
                            scale=opt$scale,
                            cluster_ids=opt$clusterids,
                            metadata_file=opt$metadata,
                            marker_table=x,
                            n_markers=24,
                            cells_use=NULL,
                            row_names_gp=8,
                        #   priority="min_log2FC",
                            sub_group=opt$group)

drawHeatmap <- function() { draw(mch) }

save_plots(paste(outprefix,"heatmap", sep="."),
           plot_fn=drawHeatmap,
           width = 7,
           to_pdf = FALSE, #opt$pdf,
           height = 2)

message("making the expression dotplots")

gp <- expressionPlots(
        loom=opt$loom,
        matrix_loc=opt$data_matrix_loc,
        gene_id_loc=opt$gene_id_loc,
        cluster_ids=opt$clusterids,
        features=x$gene_name,
        rdims=rdims,
        x=opt$rdim1,
        y=opt$rdim2,
        ncol = 4, pch = ".", 
        point_size = 1,
        max_quantile = 0.9)

gp <- gp + theme(panel.spacing.x=unit(0, "lines") ,
                 panel.spacing.y=unit(0, "lines"))

gp <- gp +theme(strip.text = element_text(size = 6, margin = margin()))

## save the plots
save_ggplots(paste(outprefix,"rdims",sep="."),
             gp,
             width=7,
             height=9,
             to_pdf=opt$pdf)

message("making the violin plots")
# make the violin plots
gg_grob <- plotHorizontalViolins(
    loom=opt$loom,
    matrix_loc=opt$data_matrix_loc,
    barcode_id_loc=opt$barcode_id_loc, #"col_attrs/barcode_id",
    gene_id_loc=opt$gene_id_loc, # "row_attrs/gene_name",
    cluster_ids=opt$clusterids,
    genes=x$gene_name,
    clusters=NULL,
    title=NULL,
    ncol=12,
    xlab="normalised expression level",
    group=NULL,
    colors=NULL,
    alpha=1,
    plot=FALSE)

do_plot <- function() { plot(gg_grob)}

do_plot()

cids <- read.table(opt$clusterids, header=T, sep="\t")

nclusters <- length(unique(cids$cluster_id))
height = min(nclusters/20 * 5,10)

## save the plots
save_plots(paste(outprefix,"violins",sep="."),
             plot_fn=do_plot,
             width=7,
             height=height,
             to_pdf=opt$pdf)

#message("making the split dotplot")

# # dot plots.
# s <- subset(s, cells=names(Idents(s)[Idents(s)==opt$cluster]))

# # sort out colors.
# if(!is.null(opt$group))
# {
#     fvar <- opt$group
#     ngroups <- length(unique(s[[fvar]][[fvar]]))
# } else {
#     fvar <- NULL
#     ngroups <- 1
# }

# cm_palette <- colormap(colormap = colormaps$portland,
#                        nshade = ngroups, alpha=0.6)

# message("building the plot")

# gp <- gp + scale_fill_manual(values=cm_palette)

# gp <- DotPlot(s,
#           features=x$gene_id,
#           cols=cm_palette,
#           split.by=fvar)

# gp <- gp + ylab(fvar) + xlab("gene")
# gp <- gp + theme_bw(base_size=8)
# gp <- gp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# message("saving the plot")
# save_ggplots(paste(outprefix,"dotplot",sep="."),
#              gp,
#              width=7,
#              height=3,
#              to_pdf=opt$pdf)


} else { message("no markers retained after filtering, skipping.")}

message("completed")
