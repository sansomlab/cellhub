## Title ----
##
## Visualisation of factors on reduced dimension (2D) plots
##
## Description ----
##
## Visualise the single-cells by their coordinates on reduced dimensions
## Cells are colored by cluster, and by the barcode fields
## which typically represent the sequence batch, sample
## and aggregation id.
##
## Details ----
##
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(ggplot2),
  require(reshape2),
  require(optparse),
  require(cellhub)
)

# Options ----

option_list <- list(
    make_option(c("--table"), default="none",
      help="A table containing the reduced coordinates and phenotype information"),
    make_option(c("--metadata"), default="none",
        help="A table containing phenotype information"),
    make_option(c("--method"),default="UMAP",
      help="Normally the type of dimension reduction"),
    make_option(c("--rdim1"), default="UMAP_1",
      help="The name of the column for reduced dimension one"),
    make_option(c("--rdim2"), default="UMAP_2",
      help="The name of the column for reduced dimension two"),
    make_option(c("--shapefactor"),default="none",
      help="A column in the cell metadata to use for deterimining the shape of points on the UMAP"),
    make_option(c("--colorfactors"),default="none",
      help="Column(s) in the cell metadata to use for deterimining the color of points on the UMAP. One plot will be made per color factor."),
    make_option(c("--analysisname", default=NULL,
      help="the name of the analysis being prented on the umap coordinates")),
    make_option(c("--plotdirvar"), default="UMAPDir",
      help="latex var holding location of plots"),
    make_option(c("--pointsize"), default=0.5,
      help="The point size for the UMAP plots"),
    make_option(c("--pointpch"), default="16",
      help="The point pch for the UMAP plots (if no shapefactor)"),
    make_option(c("--pointalpha"), default=0.8,
      help="The alpha setting for the points on the UMAP plots"),
    make_option(c("--pdf"), default=FALSE,
      help="Produce pdf versions of the plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$pointpch != ".") { opt$pointpch <- as.numeric(opt$pointpch) }

cat("Running with options:\n")
print(opt)


##
message("Reading in rdims table")
plot_data <- read.table(opt$table, sep="\t", header=TRUE)
rownames(plot_data) <- plot_data$barcode_id

if(opt$metadata!="none")
{
    meta_data <- read.table(opt$metadata, sep="\t", header=TRUE)
    mdcols <- colnames(meta_data)
    
    if("barcode_id" %in% mdcols)
    {
      rownames(meta_data) <- meta_data$barcode_id
      meta_data$barcode_id <- NULL
    } else if("barcode" %in% mdcols & "library_id" %in% mdcols)
    {
      rownames(meta_data) <- paste(gsub("-1$","",meta_data$barcode),
                                   meta_data$library_id, sep="-")
      
      meta_data$barcode <- NULL
      meta_data$library_id <- NULL
    } else {
      stop('metadata must have "barcode_id" column or "barcode" + "library_id" columns')
    }
    
    meta_cols <- colnames(meta_data)
    
    if(length(intersect(rownames(plot_data), rownames(meta_data))) < length(rownames(plot_data)))
    {
        stop("Not all cell barcodes are present in the given metadata table")
    } else {

        meta_data <- meta_data[rownames(plot_data),]
        plot_cols <- colnames(plot_data)
        plot_data <- cbind(plot_data, meta_data)
        colnames(plot_data) <- c(plot_cols, meta_cols)
    }
}



if ("cluster" %in% colnames(plot_data)){
    cluster_col = "cluster"
} else if ("cluster_id" %in% colnames(plot_data))
{
    cluster_col = "cluster_id"
} else
{
    cluster_col = "__none__"
}


if(!cluster_col=="__none__")
{
    plot_data[[cluster_col]] <- factor(plot_data[[cluster_col]],
                                       levels=sort(as.numeric(unique(plot_data[[cluster_col]]))))
}

color_vars <- strsplit(opt$colorfactors,",")[[1]]
print(color_vars)
tex = ""

message("Making UMAP plots colored by each of the factor variables")
## Make one whole-page UMAP plot per color variable
for(color_var in color_vars)
{
    print(paste("Making",color_var,"rdims plot"))

    if(color_var==cluster_col)
    {
                clust_levels = unique(plot_data[[cluster_col]])

                message("computing cluster centers")
                centers = data.frame(row.names=clust_levels, "cluster"=clust_levels,
                                     "x"=1,
                                     "y"=1)

                for(clust in clust_levels)
                {
                    centers[clust,"x"] = median(plot_data[plot_data[[cluster_col]]==clust,opt$rdim1])
                    centers[clust,"y"] = median(plot_data[plot_data[[cluster_col]]==clust,opt$rdim2])
                }
                message("computed cluster centers")

        }


    ## If a variable comprises only integers, coerce it to a character vector

    numeric = FALSE
    if(is.numeric(plot_data[[color_var]]))
    {
        if((all(plot_data[[color_var]] == round(plot_data[[color_var]]))
           & length(unique(plot_data[[color_var]])) < 50) || color_var == cluster_col)
        {
            plot_data[[color_var]] <- as.character(plot_data[[color_var]])
        }
        else
        {
            numeric = TRUE
        }

    }
    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {
        gp <- ggplot(plot_data, aes_string(opt$rdim1, opt$rdim2, color=color_var))
    } else {
        gp <- ggplot(plot_data, aes_string(opt$rdim1, opt$rdim2,
                                           color=color_var, shape=opt$shapefactor))
    }

    if(numeric)
    {
        midpoint <- (max(plot_data[[color_var]]) + min(plot_data[[color_var]]))/2
        gp <- gp + scale_color_gradient2(low="black",mid="yellow",high="red",midpoint=midpoint)
    }

    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {
        gp <- gp + geom_point(size=opt$pointsize, pch=opt$pointpch,alpha=opt$pointalpha)
    } else {
        gp <- gp + geom_point(size=opt$pointsize, alpha=opt$pointalpha)
    }


    if(color_var == cluster_col)
    {
        gp <- gp + geom_text(data=centers, aes(x, y, label=cluster), color="black", size=3)
    }


    # increase the size of the points in the legend.

    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {

        gp <- gp + guides(color = guide_legend(override.aes = list(size=5, pch=16, alpha=1)))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=5, pch=16, alpha=1)))


    } else {

        gp <- gp + guides(color = guide_legend(override.aes = list(size=5, alpha=1)))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=5, alpha=1)))


    }


    ## write out a separate legend if we have more than 20 levels.
    gp <- gp + theme_minimal()

    if(is.null(opt$analysisname))
    {
    plotfilename = paste(opt$method, color_var, sep=".")
    texCaption <- paste(opt$method,"plot colored by",color_var)
    } else {
    plotfilename = paste(opt$method, opt$analysisname, color_var, sep=".")
    texCaption <- paste(opt$method,"plot showing", opt$analysisname, "colored by",color_var)
    }
  
  
    nlevels <- length(unique(plot_data[[color_var]]))

    if(nlevels > 20 && numeric==FALSE)
    {
        legend <- g_legend(gp)
        gp <- gp + theme(legend.position="none")

        save_ggplots(file.path(opt$outdir, plotfilename),
                     gp,
                     width=7,
                     height=7,
                     to_pdf=opt$pdf)


        legendfilename = paste(opt$method, color_var, "legend", sep=".")
        save_ggplots(file.path(opt$outdir, legendfilename),
                     legend,
                     width=7,
                     height=7,
                     to_pdf=opt$pdf)

        legendCaption <- paste(texCaption, "plot legend")

        tex <- paste(tex,
                     getSubsectionTex(texCaption),
                     getFigureTex(plotfilename,texCaption,
                                  plot_dir_var=opt$plotdirvar),
                     getFigureTex(legendfilename,legendCaption,
                                  plot_dir_var=opt$plotdirvar),
                     sep="\n")

    } else {
        save_ggplots(file.path(opt$outdir, plotfilename),
                     gp,
                     width=7,
                     height=5,
                     to_pdf=opt$pdf)

        tex <- paste(tex,
                     getSubsectionTex(texCaption),
                     getFigureTex(plotfilename,texCaption,
                                  plot_dir_var=opt$plotdirvar),
                     sep="\n")
    }



}

print("Writing out latex snippet")
## write out latex snippet

tex_file <- file.path(opt$outdir,
                      paste(opt$method,
                            "tex",
                            sep="."))

writeTex(tex_file, tex)
