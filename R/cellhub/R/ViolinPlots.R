#' Get data for a set of violin plots
#' @param loom The location of the Loom file
#' @param matrix_loc The location of the matrix in the Loom file
#' @param barcode_id_loc The location of the barcodes in the Loom file
#' @param gene_id_loc The location of the gene_ids in the Loom file
#' @param scale Should the Loom data be row-scaled
#' @param cluster_ids The location of a tsv file mapping barcode_id to cluster_id
#' @param genes The set of genes to retrieve data for (row names of the seurat objects' "data" slot).
#' @param metadata A vector of meta.data columns to retrieve
#' @param clusters A vector of clusters to retrieve data for, if NULL all are included.
#'
#' @export
#'
getViolinData <- function(loom="loom.data",
                          matrix_loc="layers/log1p",
                          barcode_id_loc="col_attrs/barcode_id",
                          gene_id_loc="row_attrs/gene_name",
                          cluster_ids="cluster_ids.tsv.gz",
                          metadata_file="metadata.tsv.gz",
                          genes,
                          metadata=NULL,
                          clusters=NULL)
{
  require(reshape2)

  
  data <- getLoomData(loom_path = loom,
                   matrix_loc = matrix_loc,
                   genes = genes,
                   cells = NULL,
                   gene_id_loc = gene_id_loc,
                   barcode_id_loc = barcode_id_loc
                  )

  # this will convert hypens to dots in the barcodes.
  data <- data.frame(data)

  data$gene_name <- as.vector(rownames(data))

  message("reading in the cluster ids")
  cids <- read.table(cluster_ids, header=T, sep="\t")
  
  # see comment above re hypens.
  rownames(cids) <- gsub("-",".",cids$barcode_id)

  message("melting the data")
  ggData <- melt(data,id.vars=c("gene_name"))

  ggData$cluster <- as.numeric(cids[ggData$variable,"cluster_id"])        

  if(!is.null(clusters))
  {
    ggData <- ggData[ggData$cluster %in% clusters,]
  }

  if(length(metadata)>0)
  {
    md <- read.table(metadata_file, header=T, sep="\t")
    rownames(md) <- md$barcode_id
    ggData[,metadata] <- md[ggData$variable, metadata]
  }
  ggData
}


#' Remove spacer grobs from a facetted ggplot plot.
#' @param gp The ggplot object containing the violin plots
#' @param ncol The number of columns (as was passed to facet wrap)
#' @param nmissing  The number of facet panels missing
#'
#' @export
#'
removeGrobs <- function(gpGrob, ncol, nmissing)
{

  to_remove <- c()
  for(i in (1+(ncol-nmissing)):ncol)
  {
    to_remove <- c(to_remove, paste("panel",i,"1",sep="-"),
                   paste("strip-t",i,"1",sep="-"))
  }

  # get the grobs that must be removed
  rm_grobs <- gpGrob$layout$name %in% to_remove

  # remove grobs
  gpGrob$grobs[rm_grobs] <- NULL
  gpGrob$layout <- gpGrob$layout[!rm_grobs, ]

  gpGrob
}

#' Function to draw a block of horizontal violin plots
#'
#' The number of columns will be respected even if
#' there are fewer genes than columns
#'
#' @param ggData melted dataframe ready for plotting
#' @param title A title for this set of genes
#' @param ncol  Number of columns, will be passed to facet_wrap
#' @param group A column on which to group plots by
#' @param colors A vector of (fill) colors, one per cluster.
#' @param alpha Alpha value for the fill color
#' @param clusters Optional list of clusters to plot
#' @param cluster_labels Optional named vector containing cluster labels
#' @param xlab The label for the x axis.
#'
#' @export
#'
makeViolins <- function(ggData, 
                        title=NULL, 
                        ncol=8, 
                        group=NULL,
                        colors=NULL,
                        alpha=1,
                        clusters=NULL, cluster_labels=NULL,
                        xlab=NULL)
{
  require(ggplot2)
  require(ggstance)
  theme_set(theme_classic(base_size = 8))
  
  message("number of genes:")
  print(length(unique(ggData$gene_name)))

  ggData$gene_name <- factor(ggData$gene_name, levels=unique(ggData$gene_name))

  if(is.null(clusters))
  {
    cluster_levels <- unique(ggData$cluster)
    cluster_levels <- cluster_levels[rev(order(as.numeric(cluster_levels)))]
  } else {
    cluster_levels <- unique(clusters)
  }

  nl  <- length(levels(ggData$gene_name))
  nrow <- ceiling(nl/ncol)

  if(nl < ncol) { to_add <- (ncol - nl) } else { to_add <- 0 }

  message("preparing the data")
  if(!is.null(cluster_labels))
  {
    ggData$cluster <- factor(cluster_labels[as.character(ggData$cluster)],
                             levels=rev(as.character(cluster_labels)))
  } else {
    
    ggData$cluster <- factor(ggData$cluster,
                             levels=cluster_levels)
  }
  
  if(to_add > 0)
  {
    ggData$gene_name <- factor(ggData$gene_name,
                          levels = c(levels(ggData$gene_name),
                                     paste0(".",c(1:to_add))))
  }

  ticks <- function() {function(limits) c(round(min(limits)),max(1,floor(max(limits))))}

  message("drawing the plot")
  print(length(unique(ggData$gene_name)))
  
  if(is.null(group))
  {
    gp <- ggplot(ggData, aes(value, cluster, fill=cluster))
  } else {
    ggData[[group]] <- factor(ggData[[group]])
    gp <- ggplot(ggData, aes_string("value", "cluster", fill=group))
  }
  gp <- gp + geom_violinh(scale = "width", trim = TRUE, alpha=alpha)
  gp <- gp + facet_wrap(~gene_name, scales="free_x", ncol=ncol, drop=F)
  if(!is.null(title)) {
    gp <- gp + ggtitle(title)
  }
  if(!is.null(colors)) {
    gp <- gp + scale_fill_manual(values = colors)
  }
  gp <- gp + theme_classic()
  gp <- gp + scale_x_continuous(breaks = ticks())
  gp <- gp + xlab(xlab)
  gp <- gp + theme(legend.position = "none",
                   strip.background = element_rect(fill="grey94",color="grey94"),
                   strip.text = element_text(size=8),
                   panel.spacing.x=unit(2, "mm"),
                   axis.line.y = element_blank())

  ggData$gene_name <- droplevels(ggData$gene_name)
  
  message("adding geom_segment...")  
  gp <- gp + geom_segment(data=ggData, x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)

  message("getting the grob...") 
  gpGrob <- ggplotGrob(gp)

  if(nl < ncol) {
    message("removing grobs...") 
    gpGrob <- removeGrobs(gpGrob, ncol, ncol - nl)
  }
  message("makeViolins done.")
  gpGrob
}

#' Function to draw grob on a new page
#'
#' @param ggGrob A grob.
#'
#' @export
#'
plotGrob <- function(ggGrob)
{
  require(grid)
  grid.newpage()
  grid.draw(ggGrob)
}


#' Make a plot (or return a grob) of a set of horizontal violin plots from a seurat object.
#' @param loom The location of the Loom file
#' @param matrix_loc The location of the matrix in the Loom file
#' @param barcode_id_loc The location of the barcodes in the Loom file
#' @param gene_id_id_loc The location of the gene_ids in the Loom file
#' @param scale Should the Loom data be row-scaled
#' @param cluster_ids The location of a tsv file mapping barcode_id to cluster_id
#' @param genes The list of gene_ids to plot (should correspond to row names in the seurat object).
#' @param clusters The clusters to show - if NULL, plots are made for all clusters.
#' @param title A title for this set of genes.
#' @param ncol  Number of columns, will be passed to facet_wrap.
#' @param xlab The label for the x axis.
#' @param group A column in the seurat "meta.data" slot on which to group plots by.
#' @param colors A vector of (fill) colors, one per cluster.
#' @param cluster_labels A named vector of cluster labels
#' @param alpha Alpha level for the fill color.
#' @param plot If set to FALSE the grob is returned instead.
#'
#' @export
#'
plotHorizontalViolins <- function(loom="loom.data",
                                  matrix_loc="layers/log1p",
                                  barcode_id_loc="col_attrs/barcode_id",
                                  gene_id_loc="row_attrs/gene_name",
                                  cluster_ids="cluster_ids.tsv.gz",
                                  genes=NULL,
                                  clusters=NULL,
                                  cluster_labels=NULL,
                                  title=NULL,
                                  ncol=8,
                                  xlab="normalised expression level",
                                  group=NULL,
                                  colors=NULL,
                                  alpha=1,
                                  plot=TRUE)

{
  message("getting the data")
  ggData <- getViolinData(loom=loom,
                          matrix_loc=matrix_loc,
                          barcode_id_loc=barcode_id_loc,
                          gene_id_loc=gene_id_loc,
                          cluster_ids=cluster_ids,
                          genes=genes,
                          clusters=clusters,
                          metadata=group)

  if(!is.null(clusters) & !is.null(colors))
  {
   names(colors) <- clusters
  }
  


  message("making the violins")
  ggGrob <- makeViolins(ggData,
                        title=title,
                        ncol=ncol,
                        xlab=xlab,
                        clusters=clusters,
                        cluster_labels=cluster_labels,
                        group=group,
                        colors=colors,
                        alpha=alpha)

  
  if(plot)
  {
    message("plotting grob")
    plotGrob(ggGrob)
  } else {
    message("returning grob")
    ggGrob
  }
}
