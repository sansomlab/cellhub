require(Matrix)
require(Seurat)
require(SeuratDisk)
require(dplyr)
require(DropletUtils)

path_bender <- "/gfs/work/lkoneva/20_vincent_porcine_nuclei/analysis/imac_cellbender/output/"
samples = list.dirs(path = path_bender, recursive = F, full.names = FALSE)
samples
############## The cellbarcode are reformatted to the "UMI-LIBRARY_ID" synta.

path_mtx2 <- "/gfs/work/lkoneva/20_vincent_porcine_nuclei/analysis/imac_cellbender/mtx_umi_library_id/"

for(smp in samples){
  # load data from the filtered h5 file
  data.file <- paste0(path_bender,  smp, "/", smp, "_cellbender_filtered.h5")
  data.file
  p <- Read10X_h5(filename = data.file, use.names = FALSE)
  p
  
  # create Seurat object
  obj <- CreateSeuratObject(counts = p)
  obj
  
  head(row.names(obj))
  head(colnames(obj))
  
  feat <- read.delim("/gfs/work/lkoneva/20_vincent_porcine_nuclei/analysis/ranger/cellhub/api/cellranger.multi/GEX/filtered/sample_01/mtx/features.tsv.gz", 
                     header=FALSE)
  head(feat)
  
  # check features order
  all(rownames(obj) == feat$V1)
  identical(rownames(obj), feat$V1)
  
  # replace -1 with -smp
  head(colnames(obj))
  cell.ids <- gsub("1", smp, colnames(obj))
  head(cell.ids)
  
  write10xCounts(path = paste0(path_mtx2, smp, "/"),
                 obj@assays$RNA@counts,
                 gene.id=row.names(obj),
                 gene.symbol=feat$V2,
                 barcodes=cell.ids,
                 version="3",
                 overwrite = TRUE)
} 


###### one sample
# load data from the filtered h5 file
smp = "sample_01"
data.file <- paste0(path_cellbender, smp, "/", smp, "_cellbender_filtered.h5")
data.file
p <- Read10X_h5(filename = data.file, use.names = FALSE)
p

# create Seurat object
obj <- CreateSeuratObject(counts = p)
obj

head(row.names(obj))
head(colnames(obj))

feat <- read.delim("/gfs/work/lkoneva/20_vincent_porcine_nuclei/analysis/ranger/cellhub/api/cellranger.multi/GEX/filtered/sample_01/mtx/features.tsv.gz", 
                   header=FALSE)
head(feat)

# check features order
all(rownames(obj) == feat$V1)
identical(rownames(obj), feat$V1)

# write10xCounts(path = paste0(path_cellbender, "matrices/", smp, "/"),
#                obj@assays$RNA@counts,
#                gene.id=row.names(obj),
#                gene.symbol=feat$V2,
#                barcodes=colnames(obj),
#                version="3",
#                overwrite = TRUE)


