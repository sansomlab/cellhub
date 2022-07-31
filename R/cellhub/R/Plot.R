## function to return annotated MA plot as ggplot object
#' Make an MA-style plot
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 ratio)
#' @param a_col Column containing A value (log2 mean expression)
#' @param p_col Column containing p-values.
#' @param top_col Column indicating top rows to highlight.
#' @param label_col Column containing labels for highlighted gene_ids.
plotMA <- function(data, xlab="xlab", ylab="ylab",
                   m_col="M", a_col="A", p_col="p.adj",
                   top_col="top", label_col="gene_name",
                   p_threshold=0.05)
{

    data <- fix_p_values(data, p_col=p_col)
    data <- categoriseGenes(data,  m_col=m_col, p_col=p_col,
                            p_threshold=p_threshold)

    data <- data.frame(data)

    top <- data[data[[top_col]] == TRUE & data[[p_col]] < p_threshold,]

    ns <- data[data[[p_col]]>p_threshold,]
    sig <- data[data[[p_col]]<=p_threshold,]

    gp <- ggplot(data, aes(get(a_col), get(m_col)))
    
    gp <- gp + geom_point(data=ns, alpha=0.6, size=0.75, color="grey")
    gp <- gp + geom_point(data=sig, aes(color=-log10(get(p_col))), size=0.75 )

    gp <- gp + scale_color_viridis_c(option="rocket", direction=-1)
    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col), color="black",
                               min.segment.length=0, size=3)
    gp <- gp + geom_hline(yintercept=c(-1,1), linetype="dashed", color="grey")
    gp <- gp + geom_hline(yintercept=0, linetype="dashed")

    ltitle <- paste0("-log10\n(",p_col,")")
    gp <- gp + xlab(xlab) + ylab(ylab) + labs(col = ltitle) 

    gp
}

#' Make an volcano-style plot
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 ratio)
#' @param a_col Column containing A value (log2 mean expression)
#' @param p_col Column containing p-values.
#' @param top_col Column indicating top rows to highlight.
#' @param label_col Column containing labels for highlighted genes.
plotVolcano <- function(data, xlab="xlab", ylab="ylab",
                        m_col="M", a_col="A", p_col="p.adj",
                        top_col="top", label_col="gene_name",
                        p_threshold=0.05)
{

    data <- fix_p_values(data, p_col=p_col)
    data <- categoriseGenes(data,  m_col=m_col, p_col=p_col)
    

    data <- data.frame(data)
    top <- data[data[[top_col]] == TRUE & data[[p_col]] < 0.05,]

    ns <- data[data[[p_col]]>p_threshold,]
    sig <- data[data[[p_col]]<=p_threshold,]

    gp <- ggplot(data, aes(get(m_col), -log10(get(p_col))))

    gp <- gp + geom_point(data=ns, alpha=0.6, size=0.75, color="grey")
    gp <- gp + geom_point(data=sig, aes(color=-log10(get(p_col))), size=0.75 )
   
    gp <- gp + scale_color_viridis_c(option="rocket", direction=-1)
    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col), color="black",
                               min.segment.length=0,size=3)
    gp <- gp + geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey")
    gp <- gp + geom_vline(xintercept=c(0), linetype="dashed", color="black")
    gp <- gp + geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey")

    ltitle <- paste0("-log10\n(",p_col,")")
    gp <- gp + xlab(xlab) + ylab(ylab) + labs(col = ltitle) 

    gp
}

#' Make a plot of frequency versus expression level
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 exprs ratio)
#' @param freq_col Column containing L value (log2 freq ratio)
#' @param p_col Column containing p-values.
#' @param label_col Column containing labels for highlighted genes.
plotFvE <- function(data, p_col="p.adj", label_col="gene_name",
                    m_col="M",freq_col="L",
                    xlab="change in expression (log2)",
                    ylab="change in frequency (log2)",
                    p_threshold=0.05)
{
    data <- fix_p_values(data, p_col=p_col)
    data <- categoriseGenes(data, 
                            m_col=m_col, p_col=p_col,
                            use_fc=FALSE, ngenes=12)

    data$scolor <- "1"

    data$scolor[data$sig == TRUE] <- "2"
    data$scolor[data$top == TRUE] <- "3"

    data <- data.frame(data)
    top <- data[data$top == TRUE & data$sig,]

   
    ns <- data[data[[p_col]]>p_threshold,]
    sig <- data[data[[p_col]]<=p_threshold,]

    gp <- ggplot(data, aes(get(freq_col), get(m_col)))
    gp <- gp + geom_point(data=ns, alpha=0.6, size=0.75, color="grey")
    gp <- gp + geom_point(data=sig, aes(color=-log10(get(p_col))), size=0.75 )


    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col),
                               color="black", min.segment.length=0, size=3)

    gp <- gp + geom_hline(yintercept=0, linetype="dashed", color="grey")
    gp <- gp + geom_vline(xintercept=0, linetype="dashed", color="grey")
#   
    gp <- gp + scale_color_viridis_c(option="rocket", direction=-1)

    ltitle <- paste0("-log10\n(",p_col,")")
    gp <- gp + xlab(xlab) + ylab(ylab) + labs(col = ltitle) 

    gp
}


#' Function to make sets of violin plots.
#' @param data The data
#' @param seurat_object The seurat object
#' @param cluster_ids The cluster ids
#' @param type Either "positive" or "negative"
#' @param group.by A factor to group by
#' @param ident.include Identity to include
#' @param ncol Number of columns in the figure
#' @param m_col Column containing M value (log2 exprs ratio)
#' @param use.minfc Use minimum foldchange
#' @param minfc_col Column containing the minimum foldchange
#' @param maxfc_col Column containing the maximum foldchage
#' @param avgfc_col Column containing the average foldchage
#' @param p_col The column containing the p-value
#' @param id_col The column containing unique gene identifiers
plotViolins <- function(data=NULL, 
                        loom="loom.data",
                        matrix_loc="matrix",
                        barcode_id_loc="col_attrs/barcode_id",
                        gene_id_loc="row_attrs/gene_ids",
                        scale=FALSE,
                        cluster_ids="cluster_ids.tsv.gz",
                        type="positive",
                        group.by=NULL, ident.include=NULL,
                        vncol=4, vnrow=3,
                        use.minfc=FALSE,
                        minfc_col="min_log2FC", maxfc_col="max_logFC",
                        avgfc_col="log2FC",
                        m_col="log2FC", p_col="p.adj",
                        pt_size=0.1,
                        id_col="gene_name")
{

    ## ensure the gene ids are unique.
    data[[id_col]] <- make.unique(data[[id_col]])

    message("dimensions of data passed to plotViolins()")
    print(dim(data))

    if(type=="positive")
    {
        tmp <- data[data[[m_col]] > 0,]
    } else if (type == "negative")
    {
        tmp <- data[data[[m_col]] < 0,]
    } else {
        stop("type argument to plotViolins() not recognised")
    }

    n_show = vncol * vnrow

    ## make first sub-figure (4 columns, 3 rows)
    tmp <- tmp[order(tmp[[p_col]]),]

    genes_a <- tmp[[id_col]][min(1,length(tmp[[id_col]])):min(n_show,length(tmp[[id_col]]))]
    n_a <- length(genes_a)

    message("number of genes for first panel: ", n_a)

    nrow_a <- ceiling(n_a / vncol)

    FontSize(
        x.text = 4,
        y.text = 6,
        x.title = 0,
    #    y.title = NULL,
        main = 10
    #    ...
        )

    if(n_a > 0)
    {
        message("Calling VlnPlot for top genes by p-value")
        gpa <- VlnPlot(object=seurat_object,
                       features=genes_a,
                       pt.size = pt_size,
                       #size.x.use=6,size.y.use=6,size.title.use=10, point.size.use=0.1,
                       ncol=vncol,
                       log=TRUE,
                       # do.return=TRUE,
                       group.by=group.by,
                       idents=ident.include)

        gpa_exists <- TRUE

    } else {
        gpa_exists <- FALSE
        gpa <- FALSE
        nrow_a <- 0
    }

    ## if n_a is less than n_show, there are no more genes to show.
    if(n_a == n_show)
    {

        ## make second sub-figure
        tmp <- tmp[!tmp[[id_col]] %in% genes_a,]

        ## order by largest fold change
        ## tmp <- tmp[rev(order(abs(tmp$log2FC))),]

        ## subset to positive or negative markers
        if(type == "positive")
        {
            ## order the positive markers by fold change
            if(use.minfc)
            {
                tmp <- tmp[tmp[[minfc_col]] > 0,]
                tmp <- tmp[rev(order(tmp[[minfc_col]])),]
            } else {
                tmp <- tmp[rev(order(tmp[[avgfc_col]])),]
            }
        } else {

            ## order the negative markers by fold change
            if(use.minfc)
            {
                tmp <- tmp[tmp[[maxfc_col]] < 0,]
                tmp <- tmp[order(tmp[[maxfc_col]]),]
            } else {
                tmp <- tmp[order(tmp[[avgfc_col]]),]
            }
        }

        genes_b <- tmp[[id_col]][min(1,length(tmp[[id_col]])):min(n_show,length(tmp[[id_col]]))]
        n_b <- length(genes_b)
        nrow_b <- ceiling(n_b/vncol)

        message("number of genes for second panel: ", n_b)

        if(n_b > 0)
        {

            message("Calling VlnPlot for top genes by fold change")
            gpb <- VlnPlot(object=seurat_object,
                           features=genes_b,
                           pt.size = pt_size,
                           #size.x.use=6,size.y.use=6,size.title.use=10, point.size.use=0.1,
                           ncol=vncol,
                           log=TRUE, #do.return=TRUE,
                           group.by=group.by,
                           idents=ident.include)

            gpb_exists <- TRUE

        } else {
            gpb <- FALSE
            gpb_exists <- FALSE
            nrow_b <- 0
        }
    } else {
        gpb <- FALSE
        gpb_exists <- FALSE
        nrow_b <- 0

    }

    return(list(gpa=gpa, nrow_a=nrow_a, gpa_exists=gpa_exists,
                gpb=gpb, nrow_b=nrow_b, gpb_exists=gpb_exists))
}



#' Make a latex section containing a set of violin plots
#' @param data The data
#' @param seurat_object The seurat object
#' @param cluster_ids The cluster ids
#' @param type Either "positive" or "negative"
#' @param group.by A factor to group by
#' @param ident.include Identity to include
#' @param ncol Number of columns in the figure
#' @param use.minfc Use minimum foldchange
#' @param outdir The directory for the output
#' @param analysis_title Title for this section
#' @param fc_type The metric by which the violin plots are ordered
violinPlotSection <- function(data=NULL, 
                              loom="loom.data",
                              matrix_loc="matrix",
                              barcode_id_loc="col_attrs/barcode_id",
                              gene_id_loc="row_attrs/gene_ids",
                              cluster_ids="cluster_ids.tsv.gz",
                              #   metadata_file="metadata.tsv.gz",
                              type="positive",
                              group.by=opt$testfactor,
                              ident.include=opt$identinclude,
                              vncol=4,vnrow=3, pt_size=0.1,
                              outdir=opt$outdir,
                              analysis_title="violin plots", fc_type="fold change",
                              plot_dir_var="plotsDir",
                              to_pdf=TRUE,
                              use.minfc=FALSE)
{

    ## get the violin plots
    violin_plots <- plotViolins(data=data, 
                                loom=loom,
                                matrix_loc=matrix_loc,
                                barcode_id_loc=barcode_id_loc,
                                gene_id_loc=gene_id_loc,
                                cluster_ids=clusterids, 
                                type=type,
                                group.by=group.by, 
                                pt_size=pt_size,
                                ident.include=ident.include, 
                                vncol=vncol, 
                                vnrow=vnrow,
                                use.minfc=use.minfc)


    tex <- c()

    subsectionTitle <- getSubsectionTex(paste0("Cluster ",
                                               cluster, " violin plots: ",
                                               type,
                                               analysis_title))
    tex <- c(tex, subsectionTitle)

    ## make the figure
    tex <- c(tex, "\\begin{figure}[H]")

    if(violin_plots$gpa_exists)
    {

        ## first subfigure
        violin_fn <- gsub("type", paste0(type,".padj"), violin_fn_template)
        violin_path <- file.path(outdir, violin_fn )

        h <- max(violin_plots$nrow_a * 2, 3)

        save_ggplots(violin_path,
                     violin_plots$gpa,
                     to_pdf=to_pdf,
                     width=10,
                     height=h)

        caption <- paste0("Top ", type, analysis_title,
                          " ordered by p-value, cluster: ",cluster)

        tex <- c(tex, getSubFigureTex(violin_fn, caption, plot_dir_var=plot_dir_var, height=0.45))

        if(violin_plots$gpb_exists)
        {

            ## second subfigure
            violin_fn <- gsub("type",paste0(type,".fc"), violin_fn_template)
            violin_path <- file.path(opt$outdir, violin_fn)

            h <- max(violin_plots$nrow_b * 2, 3)

            save_ggplots(violin_path,
                         violin_plots$gpb,
                         to_pdf=to_pdf,
                         width=10,
                         height=h)

            caption <- paste0("Additional ", type, analysis_title, " ordered by ", fc_type, ", cluster: ",cluster)
            tex <- c(tex, getSubFigureTex(violin_fn, caption, plot_dir_var=plot_dir_var, height=0.45))
        }

    } else {
        tex <- c(tex, "No significant genes")
    }

    tex <- c(tex, "\\end{figure}")
    tex

}

#' Extract a legend from a ggplot
#' @param a.ggplot A ggplot obect
g_legend<-function(a.ggplot){
    tmp <- ggplot_gtable(ggplot_build(a.ggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

#' function to return colors matching ggplot2
#' see: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n Number of colours to return.
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


## save png and pdf version of base graphics plots
#' Wrapper function for saving base R plots in both pdf and png format
#' @param filepath The name of the file to write to (without .pdf or .png extension)
#' @param plot_fn The plotting function to call.
#' @param width The width of the plots.
#' @param height The height of the plots.
#' @param res The plot resolution.
#' @param to_pdf Save pdf version of plot (TRUE|FALSE).
#' @param to_png Save png version of plot (TRUE|FALSE).
save_plots <- function(filepath=NULL,
                       plot_fn=NULL,
                       width=6,
                       height=8,
                       res=300,
                       to_pdf=TRUE,
                       to_png=TRUE)
{

    ## save the pdf plot
    if(to_pdf)
    {
        cairo_pdf(paste(filepath,"pdf",sep="."),
                  width=width,
                  height=height)
        do.call(plot_fn,list())
        dev.off()

    }

    ## save the png plot
    if(to_png)
    {
        png(paste(filepath,"png",sep="."),
            width=width,
            height=height,
            units="in",
            res=res)
        do.call(plot_fn,list())
        dev.off()
    }

}


## save png and pdf versions of ggplot plots
## save png and pdf version of base graphics plots
#' Wrapper function for saving base R plots in both pdf and png format
#' @param filepath The name of the file to write to (without .pdf or .png extension)
#' @param gp The ggplot object.
#' @param width The width of the plots.
#' @param height The height of the plots.
#' @param dpi The plot resolution.
#' @param to_pdf Save pdf version of plot (TRUE|FALSE).
#' @param to_png Save png version of plot (TRUE|FALSE).
save_ggplots <- function(filepath=NULL,
                         gp=NULL,
                         width=6,
                         height=8,
                         dpi=300,
                         to_pdf=TRUE,
                         to_png=TRUE)
{
    if(to_png)
    {
        ggsave(paste(filepath,"png",sep="."),
               gp,
               device='png',
               width=width,
               height=height,
               units="in",
               type="cairo-png",
               dpi=dpi)
    }

    if(to_pdf)
    {
        ggsave(paste(filepath,"pdf",sep="."),
               gp,
               device=cairo_pdf,
               width=width,
               height=height,
               units="in",
               dpi=dpi)
    }
}

#' Plot

#' @param matrixUMI
#' @param metadata
#' @param basename

#' @return
#' @export

#' @examples
plotDownsampling <- function(matrixUMI, metadata, basename) {

    # Total UMI count per cell
    nUMIs <- Matrix::colSums(matrixUMI)

    # Collate total UMI and cell rank for each cell of each sample
    inputStats <- do.call("rbind", lapply(
        unique(metadata$sample_id),
        function(id){
            codes_id <- subset(metadata, sample_id == id, "barcode", drop=TRUE)
            nUMIs_id <- nUMIs[codes_id]
            data.frame(
                nUMIs=sort(nUMIs_id, decreasing=TRUE),
                CellRank=seq_along(nUMIs_id),
                Sample=id
            )
        }))

    # UMI vs Rank ----
    gg <- ggplot(inputStats) +
        geom_line(aes(CellRank, nUMIs, colour=Sample), size=0.25) +
        scale_x_log10(
            limits=c(1, max(inputStats$CellRank))
        ) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs))))
        ) +
        annotation_logticks() +
        labs(y="UMI count", x="Cell rank") +
        theme_bw() +
        theme(
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank()
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_ranked.pdf", basename)), gg,
        width=7, height=5
    )

    # UMI violion plot ----
    gg <- ggplot(inputStats) +
        geom_violin(
            aes(Sample, nUMIs, colour=Sample),
            draw_quantiles=c(0.5), size=0.25) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs))))
        ) +
        annotation_logticks(sides="l") +
        labs(y="UMI count", x="Sample") +
        theme_bw() +
        theme(
            # axis.text=element_text(size=rel(0.5)),
            # legend.title=element_text(size=rel(0.75)),
            # legend.text=element_text(size=rel(0.75)),
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(angle=90)
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_violin.pdf", basename)), gg,
        width=7, height=5
    )

    return(TRUE)
}


#' A function to draw a heatmap of top cluster marker genes with
#' subgroup labels. The function uses the "matrx" slot of the loom by default to make the heatmap.
#' @param loom_path A path to the loom file containing the data to be plotted
#' @param matrix_loc The location of the matrix in the loom file
#' @param barcode_id_loc The location of the barcodes
#' @param gene_id_loc The location of the genes
#' @param cluster_ids The location of a gzipped tsv file with mapping of cluster_id to barcode_id
#' @param marker_table A dataframe containing the marker information. Must contain "cluster", "gene" and "log2FC" columns
#' @param n_markers The number of markers to plot
#' @param cells_use The names of the cells to use. If NULL all cells will be used
#' @param row_names_gp The font size for the gene names
#' @param sub_group If given, the name of a variable in the metadata of the seurat object
#' @param slot scale.data|data. If data, the rows will be scaled.
#' @param scale If TRUE, rows will be scaled
#' @param disp_min Disp floor, used to truncate the scaled data
#' @param disp_max Disp ceiling, used to truncate the scaled data
#' @param padj_threshold The max acceptable padj
#' @param padj_col Name of the column with the adjusted p-values
#' @param only_positive Should only the positive markers be plotted
markerComplexHeatmap <- function(loom_path,
                                 matrix_loc="matrix",
                                 barcode_id_loc="col_attrs/barcode_id",
                                 gene_name_loc="row_attrs/gene_name",
                                 scale=FALSE,
                                 cluster_ids="cluster_ids.tsv.gz",
                                 metadata_file="metadata.tsv.gz",
                                 marker_table=NULL,
                                 n_markers=20,
                                 cells_use=NULL,
                                 padj_threshold=0.1,
                                 padj_col="p.adj",
                                 only_positive=TRUE,
                                 priority=c("pval","log2FC","min_log2FC", "score"),
                                 row_names_gp=10,
                                 sub_group=NULL,
                                 disp_min=-2.5,
                                 disp_max=2.5)
{

   # 1. Get the markers that we wish to plot
   priority <- match.arg(priority)
   rank_var <- sym(priority)
   
   # only plot significant markers
   marker_table <- marker_table[marker_table[[padj_col]] < padj_threshold & !is.na(marker_table[[padj_col]]), ]
   
   # only plot positive markers
   if(only_positive)
   {
        marker_table <- marker_table[marker_table$log2FC > 0, ]
   }
   
   if(priority=="pval")
   {
       message("sorting by p val, fold change")
       if(!"score" %in% colnames(marker_table)){
           stop("score column missing")
       }
       # sort by p val and then by fold change
       # often top markers have p==0.
   top_markers <- marker_table %>%
         group_by(cluster) %>%
         arrange(!!rank_var, desc(score)) %>%
         slice_head(n = n_markers)
   } else {
       message("sorting by ", priority)
         top_markers <- marker_table %>%
            group_by(cluster) %>%
            arrange(desc(!!rank_var)) %>%
            slice_head(n = n_markers)
   }
 
    # message("top makers gene ids:")
    # print(data.frame(top_markers[,c("gene_id", "gene_name", "cluster", "pval", "score", "log2FC")]))
    # print(table(top_markers$cluster))
  
  
  # 2. Get the data slice
    message("getting the loom data")
    x <- getLoomData(loom_path = loom_path,
                     matrix_loc = matrix_loc,
                     genes = top_markers$gene_name,
                     cells = cells_use,
                     gene_id_loc = gene_name_loc,
                     barcode_id_loc = barcode_id_loc
                     )

    if(scale) {  x <- t(scale(t(x))) }
    
    message("retrieved loom data:")
    print(dim(x))
    #print(rownames(x))

    # SWITCH OFF FOR PRODUCTION.
    # x <- x[,sample(1:ncol(x),5000)]

    # 3. sanity check the dimensions of the returned data
    gene_ids_use <- rownames(x)
    cells_use <- colnames(x)
    #print(head(cells_use))


#   } else if(length(gene_ids_use) <  length(top_markers$gene_id)) {
#      message("Warning: not all identified marker gene_ids are present in the ",
#              "selected slot. Only markers present in the selected slot ",
#              "will be plotted. If using scale.data you should consider re-scaling",
#              "your object to include all gene_ids in the scale.data slot")

#     top_markers <- top_markers %>% filter(gene_id %in% gene_ids_use)
#   }
 
  if(length(gene_ids_use) == 0) {
     stop("None of the marker gene_ids are present in the scaled data...")
  }
 
  if(length(gene_ids_use) != length(top_markers$gene_name)) {
      stop("Not all of the given genes were found in the Loom file.")
  }

  message("setting min and max values")
  # x <- MinMax(x, min = disp_min, max = disp_max)
  x[x<disp_min] <- disp_min
  x[x>disp_max] <- disp_max

  message("reading cluster assignments")
  cluster_assignments <- read.table(gzfile(cluster_ids), sep="\t",
                                    header=TRUE)
  rownames(cluster_assignments) <- cluster_assignments$barcode_id
  
  clusters <- as.numeric(cluster_assignments[cells_use, "cluster_id"])

  # set up the cluster color palette
  nclust <- length(unique(clusters))
  cluster_cols <- gg_color_hue(nclust)
  names(cluster_cols) <- sort(as.numeric(as.character(unique(clusters))))

  # get the vector of per-cell cluster names
  cell_clusters <- clusters #clusters[colnames(x)]

  message("setting the row annotations")
  clusterAnnotation = rowAnnotation(df = data.frame(cluster=top_markers$cluster),
                                    col = list(cluster=cluster_cols),
                                    show_annotation_name = FALSE,
                                    show_legend=FALSE,
                                    width=unit(2,"mm"))

  if(!is.null(sub_group))
  {
      message("setting up the subgroups")
      # ingest the metadata
      metadata <- read.table(metadata_file, header=T, sep="\t")
      rownames(metadata) <- metadata$barcode_id
      metadata <- metadata[colnames(x),]
      
      if(!sub_group %in% colnames(metadata))
      {
          stop("specified sub group not found in the metadata")
      }

    # set up the subgroup colour palette
    cell_sub_groups <- metadata[cells_use, sub_group]
    sub_groups <- unique(cell_sub_groups)

    sub_group_cols <- colormap(colormap = colormaps$portland,
                               nshade = length(sub_groups),
                               alpha = 0.6)

    # because the Grey palette returns a minimum of 3 colors..
    sub_group_cols <- sub_group_cols[1:length(sub_groups)]
    names(sub_group_cols) <- sub_groups

    print(sub_group_cols)
    # get an ordering vecotor
    cell_order <- order(cell_clusters, cell_sub_groups)

    subgroupAnnotation = HeatmapAnnotation(df=data.frame(subgroup=cell_sub_groups[cell_order]),
                                           col = list(subgroup=sub_group_cols),
                                           show_annotation_name = FALSE,
                                           annotation_legend_param = list(labels_gp = gpar(fontsize = 4),
                                                                          title_gp = gpar(fontsize = 6)))
  } else {
      
    cell_order <- order(cell_clusters)
    subgroupAnnotation <- NULL
  }

  #print(head(cell_order))
  x <- x[,cell_order]
  message("x dimensions")
  print(dim(x))
  
  message("marker table dimesions")
  print(dim(top_markers))
  
  # match the Seurat colors
  exprs_cols <- colorRamp2(c(-2.5,0,2.5), c("#FF00FF", "#000000", "#FFFF00"))

    ## fudge the row font size to something sensible.
    n_rows <- nrow(x)
    # can probably show e.g. 60 rows at font size 8
    row.cex <- min(0.5, 60/n_rows)

  # draw the heatmap
  Heatmap(x,
          cluster_rows = FALSE,
          col = exprs_cols,
          row_names_gp = gpar(fontsize = row_names_gp,
                              cex = row.cex),
          column_title_gp = gpar(fontsize = 5),
          cluster_columns = FALSE,
          show_column_names = FALSE,
          row_title = NULL,
          row_split=top_markers$cluster,
          column_split = cell_clusters[cell_order],
          top_annotation = subgroupAnnotation,
          right_annotation = clusterAnnotation,
          row_gap=unit(0.5,"mm"),
          column_title_side="bottom",
          column_gap=unit(0.5, "mm"),
          use_raster=TRUE,
          raster_device="CairoPNG",
          raster_quality=4,
          show_heatmap_legend = FALSE
  )
}


library(reshape2)


#' Draw a set of expression dotplots (2D version)
#' @param seurat_object
#' @param features A vector of feature names
#' @param rdims A data frame or matrix containing the reduced dimensions
#' @param x The name of the rdims column to use for the x coordinates
#' @param y The name of the rdims column to use for the y coordinates
#' @param ncol The number of columns to use when facet wrapping
#' @param point_size The size of the points
#' @param max_quantile The expression quantile to cap the data at.
#' @importFrom reshape2 melt
#' @export
expressionPlots <- function(
    loom="loom.data",
    matrix_loc="layers/log1p",
    barcode_id_loc="col_attrs/barcode_id",
    gene_id_loc="row_attrs/gene_name",
    cluster_ids="cluster_ids.tsv.gz",
    features=NULL,
    rdims=NULL,
    x="UMAP_1",
    y="UMAP_2",
    ncol = 6,
    pch=16,
    point_size = 2.5,
    max_quantile = 0.9) {

  require(ggplot2)
  
  cells <- rdims$barcode_id


  ncol <- min(ncol, length(features))

  message("getting data")
  data <- getLoomData(loom_path = loom,
                   matrix_loc = matrix_loc,
                   genes = features,
                   cells = cells,
                   gene_id_loc = gene_id_loc,
                   barcode_id_loc = barcode_id_loc
                  )
  
  # transform to % of 90th quantile so that a common
  # colour scale can be used.
  scaled_data <- apply(data, 1, FUN = scale_to_quantile, q = max_quantile)

  fill_frame <-cbind(rdims[,c(x, y)],
                     scaled_data)

  fill_df <- melt(fill_frame, id.vars=c(x, y))

  gp <- ggplot(fill_df, aes_string(x, y, color="value"))
  gp <- gp + geom_point(size=point_size,
                        alpha=1,
                        stroke = 0,
                        shape = pch)
  gp <- gp + scale_color_gradientn(colours=c("grey","yellow","red"))
  gp <- gp + facet_wrap(~variable, ncol= ncol)
  gp <- gp + theme_minimal(base_size=6)
  gp

}


#' Draw a set of expression dotplots (3D version)
#' @param seurat_object
#' @param features A vector of feature names
#' @param rdims A data frame or matrix containing the reduced dimensions
#' @param x The name of the rdims column to use for the x coordinates
#' @param y The name of the rdims column to use for the y coordinates
#' @param z The name of the rdims column to use for the z coordinates
#' @param ncol The number of columns to use when facet wrapping
#' @param point_size The size of the points
#' @param max_quantile The expression quantile to cap the data at.
#' @importFrom reshape2 melt
#' @export
expressionPlots3D <- function(seurat_object, features,
                              rdims,
                              x="UMAP_1",
                              y="UMAP_2",
                              z="UMAP_3",
                              ncol = 6,
                              point_size=2.5,
                              theta=0,
                              phi=130,
                              draw_axes=FALSE,
                              max_quantile=0.9) {
  require(ggplot2)
  require(gg3D)

  cells <- rdims$barcode

  checkFeatures(seurat_object, features)
  checkCells(seurat_object, cells)

  ncol <- min(ncol, length(features))

  data <- GetAssayData(seurat_object, slot =  "data",
                       assay="RNA")[features, cells]

  # transform to % of 90th quantile so that a common
  # colour scale can be used.

  scaled_data <- apply(data, 1, FUN = scale_to_quantile, q = max_quantile)


  fill_frame <-cbind(rdims[,c(x, y, z)],
                     scaled_data)


  fill_df <- melt(fill_frame, id.vars=c(x, y, z))


  print("drawing the plot")

  gp <- ggplot(fill_df, aes_string(x=x, y=y, z=z, color="value"))
  gp <- gp + theme_void()

  if(draw_axes) {
  gp <- gp + axes_3D(theta=theta, phi=phi)
  }

  gp <- gp + stat_3D(theta=theta, phi=phi, geom="point", alpha=1,
                     size=point_size,
                     stroke=0, shape=16)

  gp <- gp + scale_color_gradientn(colours=c("grey","yellow","red"))
  gp <- gp + facet_wrap(~variable, ncol= ncol)

  gp
}
