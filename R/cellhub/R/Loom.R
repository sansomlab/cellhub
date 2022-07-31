## function to open a connection to a loom file
#' @param file_path the path to the loom file
#' @param mode 
#' @param matrix_loc the location of the matrix
#' @param gene_id_loc the location of the gene_ids
#' @param barcode_id_loc the location of the cell barcodes
connectLoom <- function(file_path, 
                        mode='r',
                        matrix_loc='matrix',
                        gene_id_loc="row_attrs/gene_ids",
                        barcode_id_loc="col_attrs/barcode_id")
{
  loomCon <- loomR::connect(filename = file_path, 
                     mode = mode,
                     skip.validate=TRUE)
  
  message("checking that required attributes and matrices are present")
  check <- loomCon[[barcode_id_loc]][1]
  check <- loomCon[[gene_id_loc]][1]
  check <- loomCon[[matrix_loc]][0,0]
  
  message("returning the loom connection")
  return(loomCon)
}

## function to extract the numeric indices
## of cells and genes of interest from a loom file
#' @param loomCon the loom file connection
#' @param genes the gene identifiers for which to retrieve indices
#' @param cells the cell identifiers for which to retrieve indices
#' @param gene_id_loc the location of the gene_ids
#' @param barcode_id_loc the location of the cell barcodes
getLoomIndices <- function(loomCon=NULL,
                           genes=NULL, 
                           cells=NULL,
                           gene_id_loc="row_attrs/gene_ids", 
                           gene_name_loc="row_attrs/gene_name",
                           barcode_id_loc="col_attrs/barcode_id")
{
  
   # make a hash structure for mapping barcodes to indices
  nbarcodes <- loomCon[[barcode_id_loc]]$dims
  barcode_hash <- 1:nbarcodes
  names(barcode_hash) <- loomCon[[barcode_id_loc]][barcode_hash]
  
  # make a hash structure for mapping genes to indices
  ngenes <- loomCon[[gene_id_loc]]$dims
  gene_hash <- 1:ngenes
  names(gene_hash) <- loomCon[[gene_id_loc]][gene_hash]

  if(is.null(genes))
  {
    message("Getting indicies for all genes present in loom")
    genes <- names(gene_hash)
  } else {
     ngenes_requested <- length(genes)
     message("No. genes requested:", ngenes_requested )
     # print(genes)
     genes <- genes[genes %in% names(gene_hash)]
    
     ngenes_found <- length(genes)
     message("No. genes found:", ngenes_found)
     # print(genes)
     
  
    if(ngenes_found < ngenes_requested)
    {
      warning("Not all genes requested were found in the loom.")
    } 
  }
  
  if(is.null(cells))
  {
    message("Getting indices for all cells present in loom")
    cells <- names(barcode_hash)
  } else {
    
    ncells_requested <- length(cells)
    cells <- cells[cells %in% names(barcode_hash)]
    ncells_found <- length(cells)
    message("No. cells requested:", ncells_requested )
    message("No. cells found:", ncells_found)
    if(ncells_found < ncells_requested)
    {
      warning("Not all cells requested were found in the loom.")
    }
  }
  
  cellix <- as.vector(barcode_hash[cells])
  cell_names <- names(barcode_hash[cells])
  geneix <- as.vector(gene_hash[genes])
  gene_ids <- names(gene_hash[genes])
  gene_names <- loomCon[[gene_name_loc]][geneix]
  
  return(list(cells=cellix, 
              cell_names=cell_names,
              genes=geneix,
              gene_ids=gene_ids,
              gene_names=gene_names))
}



## function to extract the numeric indices
## of cells and genes of interest from a loom file
#' @param loomCon the loom file connection
#' @param genes the gene identifiers for which to retrieve indices
#' @param cells the cell identifiers for which to retrieve indices
#' @param gene_ids the location of the gene_ids
#' @param barcodes the location of the cell barcodes
getLoomData <- function(loom_path=NULL,
                        matrix_loc="matrix",
                        gene_id_loc="row_attrs/gene_ids", 
                        gene_name_loc="row_attrs/gene_name", 
                        barcode_id_loc="col_attrs/barcode_id",
                        genes=NULL, 
                        cells=NULL)
{
  
  if(!file.exists(loom_path)) { stop("invalid loom file path")}

  # get the data slice from the loom file
  message("connecting to loom")
  loomCon = connectLoom(file_path=loom_path,
                        matrix=matrix_loc,
                        gene_id_loc=gene_id_loc,
                        barcode_id_loc=barcode_id_loc
                        )

  message("fetching barcode and gene indices")
  ix <- getLoomIndices(loomCon=loomCon,
                       genes=genes, 
                       cells=cells,
                       gene_id_loc=gene_id_loc,
                       gene_name_loc=gene_name_loc,
                       barcode_id_loc=barcode_id_loc)


  message("fetching the data from the loom")
  x <- as.matrix(loomCon[[matrix_loc]][ix$cells, 
                                       ix$genes])

  colnames(x) <- ix$gene_ids
  rownames(x) <- ix$cell_names
  
  # transpose so cells are in columns.
  x <- t(x)

  x
  
}


