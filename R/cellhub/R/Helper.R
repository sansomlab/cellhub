## Generic helper functions

#' Tidy up numbers in results tables
#' @param results_table A dataframe or matrix of results contaning numeric data.
#' @param nsignif Number of significant figures (used for numbers > 1)
#' @param nround Number of digits (used for numbers < 1)
tidyNumbers <- function(results_table, nsignif=3, nround=2)
{
  x <- as.data.frame(results_table)
  for(col in colnames(x))
  {
    if(is.numeric(x[[col]]))
    {
      a <- abs(x[[col]])< 1 & !is.na(x[[col]])
      b <- abs(x[[col]])> 1 & !is.na(x[[col]])
      x[[col]][a] <- signif(x[[col]][a],nsignif)
      x[[col]][b] <- round(x[[col]][b],nround)
    }
  }
  x
}

#' format numbers in results tables as character strings
#'
#' @param results_table A dataframe or matrix of results that contains some numeric columns.
#' @param number_fmt Sprintf format for numeric columns
sprintfResults <- function(results_table,
                           number_fmt="%0.3g")
{
    for(col in colnames(results_table))
    {
        if(is.numeric(results_table[[col]]))
        {
            xx <- results_table[[col]]

            nas <- is.na(xx)

            xx[!nas] <- sapply(xx[!nas],
                               sprintf,
                               fmt=number_fmt)

            results_table[[col]] <- xx
        }
    }

    results_table
}


#' Function to add "top" indicator column to degenes
#' @param data The data
#' @param m_col The column containing the log2 ratio
#' @param id_col A column containing unique identifiers
#' @param ngenes The number of genes to demarcate
topGenes <- function(data, m_col = "log2FC",
                     use_fc = TRUE,
                     id_col="gene", ngenes=7)
{

  tmp <- data[order(data$p.adj),id_col][1:ngenes]
  if(use_fc)
  {
    # ensure non-redunant
    data <- data[!data[[id_col]] %in% tmp,]
    tmp <- c(tmp, data[rev(order(abs(data[[m_col]]))),id_col][1:ngenes])
  }
  tmp <- unique(tmp)
  tmp
}

#' Function to add "top" and "sig" indicator columns to degenes list
#' @param data The data
#' @param m_col The column containing the log2 ratio
#' @param p_col The column containing the p-value
#' @param p_threshold The p-threshold below which genes are considered significant
#' @param m_col The column containing the log2 ratio
#' @param ngenes The number of genes to demarcate
#' @param id_col A column containing a unique identifier
categoriseGenes <- function(data,m_col="log2FC", use_fc=TRUE,
                            p_col="p.adj", p_threshold=0.05,
                            ngenes=7,
                            id_col="gene_name")
{

  tmp <- topGenes(data[data[[m_col]] > 0,],
                  m_col=m_col, use_fc=use_fc,
                  ngenes=ngenes,id_col=id_col)
  tmp2 <- topGenes(data[data[[m_col]] < 0,],
                   m_col=m_col, use_fc=use_fc,
                   ngenes=ngenes,id_col=id_col)
  data$top <- FALSE
  data$top[data[[id_col]] %in% unique(c(tmp,tmp2))] <- TRUE


  data$sig <- FALSE
  data$sig[data[[p_col]] < p_threshold] <- TRUE
  data
}

#' Function to get the significant principle components following
#' a jackstraw analysis
#' @param ncomp If set, the max. number of significant components to return
#' @param adjust_p Should p values be adjusted
#' @param p_adjust_method Method for adjusting p values
#' @param p_threshold The significance threshold
#' @param order_by_sig If true, PCs are selected based on significance
#' @param seurat_object A seurat object on which JackStraw and JackStrawPlot have been run
getSigPC <- function(seurat_object=NULL,
                     ncomp=Inf,
                     adjust_p=TRUE,
                     p_adjust_method="BH",
                     p_threshold=0.05,
                     order_by_sig=FALSE)
{
    x <- as.data.frame(seurat_object@reductions$pca@jackstraw@overall.p.values)

    ## adjust the p values
    if(adjust_p)
        {
            x$p <- p.adjust(x$Score, method=p_adjust_method)
        } else { x$p <- x$Score}

    ## subset to significant
    x <- x[x$p < p_threshold,]

    ## reorder by significance
    if(order_by_sig)
    {
        x<-x[order(x$p),]
    }

    ## return the significant components (ensuring logical order...)
    ntake <- min(nrow(x),ncomp)

    x$PC[1:ntake]

}

#' A function to check if features are present in a Seurat object
#' @param seurat_object A seurat object
#' @param features A vector of the feature names to be checked
checkFeatures <- function(seurat_object, features) {
  
  if(!all(features %in% rownames(seurat_object)))
  {
    print("these features were not found in the seurat object:")
    print(features[!features %in% rownames(seurat_object)])
    stop()
  }
  
}

#' A function to check if cells are present in a Seurat object
#' @param seurat_object A seurat object
#' @param features A vector of the feature names to be checked
checkCells <- function(seurat_object, cells) {
  
  if(!all(cells %in% Cells(seurat_object)))
  {
    missing = cells[!cells %in% Cells(seurat_object)]
    print(paste("Some of the cells were not found in the seurat object, n=",
                length(missing),
                sep=""))
    print("Examples of the missing cells:")
    print(head(missing))
    stop()
  }
  
}

#' Rescale an expression matrix as a percentage of the given quantile 
#' @param x the data
#' @param q the quantile to use as the maximum value
scale_to_quantile <- function(x,q=0.9)
{
  # cap at an upper quantile
  qv <- quantile(x[x!=0],q)
  x[x >= qv] <- qv
  
  # transform to %scale 
  x <- x/qv
  x
}

#'  sort out the p-values for plotting purposes
#'  sometimes p-values will be set to NA if not computed
#'  zeros are a problem for plotting
#' @param data The data
#' @param p_col The p value column in the data
fix_p_values <- function(data, p_col="p.adj")
   {
    data[[p_col]][is.na(data[[p_col]])] <- 1
    
    # zeros generate Infs and are excluded from the color scale.
    min_p <- data[[p_col]][data[[p_col]]>0]
    # set the zeros to a very small value.
    min_p <- min(1e-16, min_p)
    data[[p_col]][data[[p_col]]==0] <- min_p
    
    data
   }

#' This function returns a map of ensembl_id -> gene_name
#' @param annotation_file A tsv file with mapping of ensembl_id to gene_name
parseBiomartAnnotation <- function(annotation_file)
{
  anno <- read.table(annotation_file, header=T, sep="\t")
  anno <- unique(anno[,c("ensembl_id","gene_name")])
  
  missing_names <- anno$gene_name==""
  anno$gene_name[missing_names] <- anno$ensembl_id[missing_names]
  
  # make the gene names unique
  x <- make.unique(anno$gene_name)
  names(x) <- anno$ensembl_id
  x
}

#' Map ensembl_ids to gene_names with a given map
getGeneNames <- function(anno_map, ensembl_ids)
{
  if(!all(ensembl_ids %in% names(anno_map)))
  {
    warning("The map is missing some of the ensembl ids")
    message("number of ensembl ids not present in the map:")
    print(length(ensembl_ids[!ensembl_ids %in% names(anno_map)]))
    message("adding missing to map (with id as gene_name)")
    extra <- unique(ensembl_ids[!ensembl_ids %in% names(anno_map)])
    names(extra) <- extra
    anno_map <- c(anno_map, extra)
  }
  x <- as.vector(anno_map[ensembl_ids])
  x
}
