# create exclude lists of contaminant cells


# ----------
# Myeloid 
# ---------

# cluster results

cls_dir <- "~/work/synovium_atlas/cellhub/allstudies/myeloid/scxl/myeloid.seurat.dir/components.15.dir"

cl_dirs <- list.dirs(path = cls_dir, recursive = FALSE)
cl_dirs <- cl_dirs[grep("cluster.", cl_dirs)]
cl_dirs
names(cl_dirs) <- gsub(".dir", "", basename(cl_dirs))
print(cl_dirs)

clusters <- lapply(names(cl_dirs), function(cl) {
  cl_dir <- cl_dirs[[cl]]
  cl_tab <- fread(paste0(cl_dir, '/', 'scanpy.clusters.tsv.gz'))
  colnames(cl_tab)[2] <- cl
  cl_tab
})

cluster_tab <- Reduce(function(x, y) merge(x, y, by = "barcode", sort = FALSE), 
                      clusters)
head(cluster_tab)
# cluster 28 resolution 2: fibroblast contamination

exclude <- filter(cluster_tab, cluster.2 == 28) %>%
  select(barcode)
head(exclude)

write.table(exclude, "~/work/synovium_atlas/cellhub/allstudies/myeloid/exclude_barcode_list_fibro.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ----------
# Vascular cells 
# ---------

# cluster results

cls_dir <- "~/work/synovium_atlas/cellhub/allstudies/vascularcell/scxl/vascularcell.seurat.dir/components.15.dir"

cl_dirs <- list.dirs(path = cls_dir, recursive = FALSE)
cl_dirs <- cl_dirs[grep("cluster.", cl_dirs)]
cl_dirs
names(cl_dirs) <- gsub(".dir", "", basename(cl_dirs))
print(cl_dirs)

clusters <- lapply(names(cl_dirs), function(cl) {
  cl_dir <- cl_dirs[[cl]]
  cl_tab <- fread(paste0(cl_dir, '/', 'scanpy.clusters.tsv.gz'))
  colnames(cl_tab)[2] <- cl
  cl_tab
})

cluster_tab <- Reduce(function(x, y) merge(x, y, by = "barcode", sort = FALSE), 
                      clusters)
head(cluster_tab)
# cluster 15 resolution 2: fibroblast contamination

exclude <- filter(cluster_tab, cluster.2 == 15) %>%
  select(barcode)
head(exclude)
dim(exclude)

write.table(exclude, "~/work/synovium_atlas/cellhub/allstudies/vascularcell/exclude_barcode_list_fibro.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ----------
# T cells 
# ---------

# cluster results

cls_dir <- "~/work/synovium_atlas/cellhub/allstudies/vascularcell/scxl/vascularcell.seurat.dir/components.15.dir"

cl_dirs <- list.dirs(path = cls_dir, recursive = FALSE)
cl_dirs <- cl_dirs[grep("cluster.", cl_dirs)]
cl_dirs
names(cl_dirs) <- gsub(".dir", "", basename(cl_dirs))
print(cl_dirs)

clusters <- lapply(names(cl_dirs), function(cl) {
  cl_dir <- cl_dirs[[cl]]
  cl_tab <- fread(paste0(cl_dir, '/', 'scanpy.clusters.tsv.gz'))
  colnames(cl_tab)[2] <- cl
  cl_tab
})

cluster_tab <- Reduce(function(x, y) merge(x, y, by = "barcode", sort = FALSE), 
                      clusters)
head(cluster_tab)
# cluster 24 resolution 2: PRG4 fibroblast contamination
# cluster 32 resolution 3: cycling B-cell contamination

exclude <- filter(cluster_tab, cluster.2 == 15) %>%
  select(barcode)
head(exclude)
dim(exclude)

write.table(exclude, "~/work/synovium_atlas/cellhub/allstudies/vascularcell/exclude_barcode_list_fibro.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


