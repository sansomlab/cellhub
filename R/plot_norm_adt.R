# Libraries ------
stopifnot(require(optparse),
          require(futile.logger),
          require(dplyr),
          require(data.table),
          require(Matrix),
          require(Seurat),
          require(ggplot2),
          require(cowplot)
)

# Global options ------

options(stringsAsFactors = F)

# Script arguments ------

option_list <- list(
  make_option(c("--cellranger_dir"), default=".",
              help="Folder with filtered cellranger output. Must include 
              barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--library_id"), default=NULL,
              help="library or channel id"),
  make_option(c("--gex_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--adt_depth"), default=NULL,
              help="Table with the automatic prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--bcmin"), default=NULL,
              help="ADT UMI count mininum limit of the background interval"),
  make_option(c("--bcmax"), default=NULL,
              help="ADT UMI count maximum limit of the background interval"),
  make_option(c("--bfmin"), default=NULL,
              help="GEX mininum number of detected features of the background 
              interval"),
  make_option(c("--bfmax"), default=NULL,
              help="GEX maxinum number of detected features of the background 
              interval"),
  make_option(c("--ccmin"), default=NULL,
              help="ADT UMI count mininum limit of the cell-containing 
              interval"),
  make_option(c("--ccmax"), default=NULL,
              help="ADT UMI count maxinum limit of the cell-containing 
              interval"),
  make_option(c("--cfmin"), default=NULL,
              help="GEX mininum number of detected features of the 
              cell-containing interval"),
  make_option(c("--cfmax"), default=NULL,
              help="GEX maxinum number of detected features of the 
              cell-containing interval"),
  make_option(c("--numcores"), default=2,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="plot_norm_adt.log"),
  make_option(c("--outfile"), default="./report_norm_adt.pdf")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger ------

# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -
###################
flog.threshold(INFO)

# now set to append mode -------------------------------------------------------
###################
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

# Read in data -----------------------------------------------------------------
###################

flog.info("Reading data object...")
raw_adt_mtx <- Read10X(file.path(opt$cellranger_dir))
flog.info("Number of barcodes in input: %s", format(ncol(raw_adt_mtx), 
                                                    big.mark=","))

gex_depth <- fread(opt$gex_depth)
adt_depth <- fread(opt$adt_depth)

# Visualize unfiltered ADT count matrix into bacground & cell containing 
# matrices ---------------------------------------------------------------------
###################

get_den <- function(x, v="log_total_UMI") {
  ggplot(x, aes(!!rlang::sym(v))) +
    geom_density() +
    theme_classic() +
    ylab("Density")
}

get_dens <- function(x, v="log_total_UMI", g="group") {
  ggplot(x, aes(!!rlang::sym(v))) +
    geom_density(aes(fill = !!rlang::sym(g)), alpha = .5) +
    theme_classic() +
    ylab("Density") +
    theme(legend.position = "bottom")
}

flog.info("Defining background and cell containing barcodes based on \
            authomated procedure.")
  
# ADT background vs cell plots
filter(adt_depth, `group` == "background")[["log_total_UMI"]] %>%
  min -> bcmin
filter(adt_depth, `group` == "background")[["log_total_UMI"]] %>%
  max -> bcmax
  
gg_adt_den <- get_den(adt_depth, "log_total_UMI") +
  geom_vline(xintercept = c(bcmin, bcmax)) +
  ggtitle("ADT background automated")
  
sub_adt_depth <- filter(adt_depth,
                        group %in% c("cell",
                                     "background"))
  
gg_adt_dens <- get_dens(sub_adt_depth) +
  ggtitle(paste("# Back_Cell barcodes:", table(sub_adt_depth$group), 
          collapse = "_"))
  
# GEX background vs cell plots
filter(gex_depth, `group` == "background")[["ngenes"]] %>%
  min -> bfmin
  
filter(gex_depth, `group` == "background")[["ngenes"]] %>%
  max -> bfmax
  
gg_gex_den <- get_den(gex_depth, "ngenes") +
  geom_vline(xintercept = c(bfmin, bfmax)) +
  ggtitle("GEX background automated") +
  scale_x_log10() 
  
sub_gex_depth <- filter(gex_depth,
                        group %in% c("cell",
                                     "background"))
  
gg_gex_dens <- get_dens(sub_gex_depth, "ngenes") +
  facet_wrap(~group, scale = "free") +
  ggtitle(paste("# Back_Cell barcodes:", table(sub_gex_depth$group), 
          collapse = "_"))

# Intersection
gex_cell_barcodes <- filter(gex_depth, group == "cell")[["barcode_id"]]
adt_cell_barcodes <- filter(adt_depth, group == "cell")[["barcode_id"]]
cell_barcodes <- intersect(gex_cell_barcodes, adt_cell_barcodes)

flog.info("Number of cell barcodes:", length(cell_barcodes), 
          capture = TRUE)

gex_back_barcodes <- filter(gex_depth, group == "background")[["barcode_id"]]
back_barcodes <- gex_back_barcodes[!gex_back_barcodes %in% cell_barcodes]

# ADT final background vs cell plots
inter_adt_depth <- filter(adt_depth, 
                            barcode_id %in% cell_barcodes |
                            barcode_id %in% back_barcodes) %>%
    mutate('group' = ifelse(barcode_id %in% cell_barcodes, "cell",
                            ifelse(barcode_id %in% back_barcodes, "background",
                                   "???")))
  
gg_int_adt_den <- get_den(inter_adt_depth, "log_total_UMI") +
  ggtitle("ADT background vs cell automated")
  
gg_int_adt_dens <- get_dens(inter_adt_depth) +
  ggtitle(paste(table(inter_adt_depth$group), 
                collapse = "_"))
  
pdf(opt$outfile, width = 8, height = 14)
  plot_grid(plotlist = list(gg_adt_den,
                            gg_adt_dens,
                            gg_gex_den,
                            gg_gex_dens,
                            gg_int_adt_den,
                            gg_int_adt_dens), 
            ncol = 2)
dev.off()
  
flog.info("Number of background barcodes:", length(back_barcodes), 
          capture = TRUE)
  
if(opt$bcmin != "None") {
  
  flog.info("Ploting background and cell containing barcodes based on \
            user's thresholds")
  
  gex_cell_barcodes <- filter(gex_depth,
                              ngenes > opt$cfmin &
                                ngenes < opt$cfmax)[["barcode_id"]]
  
  adt_cell_barcodes <- filter(adt_depth,
                              log_total_UMI > opt$ccmin &
                                log_total_UMI < opt$ccmax)[["barcode_id"]]
  
  cell_barcodes <- intersect(gex_cell_barcodes, adt_cell_barcodes)
  
  flog.info("Number of cell barcodes:", length(cell_barcodes), 
            capture = TRUE)
  
  adt_back_barcodes <- filter(adt_depth,
                              log_total_UMI > opt$bcmin &
                                log_total_UMI < opt$bcmax)[["barcode_id"]]
  
  back_barcodes <- adt_back_barcodes[!adt_back_barcodes %in% cell_barcodes]
  
  # ADT user's defined background vs cell plots
  inter_adt_depth <- filter(adt_depth, 
                            barcode_id %in% cell_barcodes |
                              barcode_id %in% back_barcodes) %>%
    mutate('group' = ifelse(barcode_id %in% cell_barcodes, "cell",
                            ifelse(barcode_id %in% back_barcodes, "background",
                                   "???")))
  
  gg_user_adt_den <- get_den(inter_adt_depth, "log_total_UMI") +
    ggtitle("ADT user's defined background vs cells")
  
  gg_user_adt_dens <- get_dens(inter_adt_depth) +
    ggtitle(paste(table(inter_adt_depth$group), 
                  collapse = "_"))
  
  pdf(opt$outfile, width = 8, height = 18)
  plot_grid(plotlist = list(gg_adt_den,
                            gg_adt_dens,
                            gg_gex_den,
                            gg_gex_dens,
                            gg_int_adt_den,
                            gg_int_adt_dens,
                            gg_user_adt_den,
                            gg_user_adt_dens), 
            ncol = 2)
  dev.off()
  
  flog.info("ADT normalization report written.")
  
}