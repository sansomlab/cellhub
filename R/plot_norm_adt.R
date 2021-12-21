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
  make_option(c("--unfiltered_dir"), default=".",
              help="Folder with unfiltered cellranger output. Must include
              barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz"),
  make_option(c("--library_id"), default=NULL,
              help="library or channel id"),
  make_option(c("--gex_depth"), default=NULL,
              help="Table with the cellranger prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--adt_depth"), default=NULL,
              help="Table with the cellranger prediction of GEX bacground and 
              cell-containig barcodes"),
  make_option(c("--cap_thr"), default=0.02,
              help="Top fraction of most & last cell-barcodes to remove from 
              analysis. (Number of genes, Number of UMI.)"),
  make_option(c("--rm_feat"), default="None",
              help="Comma separated string with features to remove from
              analysis."),
  make_option(c("--qc_bar"), default="None",
              help="Single column .tsv.gz file with the cell-barcodes with high
              GEX-based QC metrics. (Output fetch_cells pipeline.)"),
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
raw_adt_mtx <- Read10X(file.path(opt$unfiltered_dir))
flog.info("Number of barcodes in input: %s", format(dim(raw_adt_mtx), 
                                                    big.mark=","))

if(opt$rm_feat != "None") {
  rmft <- unlist(strsplit(opt$rm_feat, ","))
  if(any(!rmft %in% rownames(raw_adt_mtx))) {
    stop("Features provided to be removed are not present.")
  }
  raw_adt_mtx <- raw_adt_mtx[!rownames(raw_adt_mtx) %in% rmft, ]
}

flog.info("Number of barcodes in input: %s", format(dim(raw_adt_mtx), 
                                                    big.mark=","))

gex_depth <- fread(opt$gex_depth)
flog.info("GEX depth:", head(gex_depth), capture = TRUE)

if(opt$qc_bar != "None") {
  qcbar <- fread(opt$qc_bar, h = F)[[1]]
  if(length(which(qcbar %in% gex_depth$barcode_id)) == 0) {
    stop("QC barcodes are not present.")
  }
  gex_depth <- filter(gex_depth, `group` == "background" | barcode_id %in% qcbar)
}

adt_depth <- fread(opt$adt_depth)
flog.info("ADT depth:", head(adt_depth), capture = TRUE)

if(opt$qc_bar != "None") {
  qcbar <- fread(opt$qc_bar, h = F)[[1]]
  if(length(which(qcbar %in% adt_depth$barcode_id)) == 0) {
    stop("QC barcodes are not present.")
  }
  adt_depth <- filter(adt_depth, `group` == "background" | barcode_id %in% qcbar)
}

# Visualize unfiltered ADT count matrix into background & cell containing 
# matrices ---------------------------------------------------------------------

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
            cellranger procedure.")
  
# ADT background vs cells. Plots with no filter --------------------------------

quantile(filter(adt_depth, `group` == "cell")[["log_total_UMI"]], 
         opt$cap_thr, na.rm = TRUE) -> ccmin
quantile(filter(adt_depth, `group` == "cell")[["log_total_UMI"]], 
         1-opt$cap_thr, na.rm = TRUE) -> ccmax

if(ccmin == 0) {
  quantile(filter(adt_depth, `group` == "cell")[["log_total_UMI"]], 
           opt$cap_thr + .02, na.rm = TRUE) -> ccmin
}

flog.info("Cell ADT limits (ccmin: cell-count min value):", ccmin, 
          capture = TRUE)
flog.info("Cell ADT limits (ccmin: cell-count max value):", ccmax, 
          capture = TRUE)

gg_adt_den <- get_den(adt_depth, "log_total_UMI") +
  geom_vline(xintercept = c(ccmin, ccmax)) +
  ggtitle("ADT background cellranger",
          subtitle = "Cell limits (96% of barcodes)")
  
sub_adt_depth <- filter(adt_depth,
                        group %in% c("cell",
                                     "background"))
  
gg_adt_dens <- get_dens(sub_adt_depth) +
  ggtitle("ADT background vs cell (cellranger)",
          subtitle = paste(table(sub_adt_depth$group),
                           collapse = "_")) +
  xlab("ADT log10(UMI)")

# GEX background vs cell plots -------------------------------------------------

quantile(filter(gex_depth, `group` == "cell")[["nfeat"]], opt$cap_thr, 
         na.rm = TRUE) -> cfmin
quantile(filter(gex_depth, `group` == "cell")[["nfeat"]], 1-opt$cap_thr,
         na.rm = TRUE) -> cfmax

if(cfmin < 200) {
  flog.info("WARNING!: Minimum number of genes in cell is less than 200", cfmin, 
            capture = TRUE)
} 

flog.info("Cell GEX limits (cfmin: cell-feature min value):", cfmin, 
          capture = TRUE)
flog.info("Cell GEX limits (cfmin: cell-feature max value):", cfmax, 
          capture = TRUE)

gg_gex_den <- get_den(gex_depth, "nfeat") +
  geom_vline(xintercept = c(cfmin, cfmax)) +
  ggtitle("GEX background cellranger",
          subtitle = paste0("Cell limits (", ceiling((1-(opt$cap_thr)*2)*100), 
                            "% of barcodes)")) +
  xlab("# genes") +
  scale_x_log10() 
  
sub_gex_depth <- filter(gex_depth,
                        group %in% c("cell",
                                     "background"))
  
gg_gex_dens <- get_dens(sub_gex_depth, "nfeat") +
  facet_wrap(~group, scale = "free") +
  ggtitle("GEX background vs cell (cellranger)",
          subtitle = paste(table(sub_gex_depth$group), 
                           collapse = " ")) +
  xlab("# genes")

# Plots trimming background ---------------------------------------------------

sub_gex_depth <- gex_depth[, c(1, 6, 7)]
colnames(sub_gex_depth)[2] <- "log_nfeat_gex"
head(sub_gex_depth)

subadt_depth <- adt_depth[, c(1, 5)]
colnames(subadt_depth)[2] <- "log_total_UMI_adt"
head(subadt_depth)

depth <- merge(sub_gex_depth, 
               subadt_depth, 
               by = "barcode_id")

head(depth)
flog.info("Dim depth:", dim(depth), capture = TRUE)
flog.info("Depth classes:", table(depth$group), capture = TRUE)

depth %>%
  filter(log_total_UMI_adt >= 0.5) -> sub_depth
  #filter(log_nfeat_gex >= 0.5) 

flog.info("Dim sub-depth:", dim(sub_depth), capture = TRUE)

back <- filter(sub_depth, group == "background")
flog.info("Barcode dimensions:", dim(back), capture = TRUE)

if(nrow(back) < 100) {
  
  flog.info("WARNING: Number of background barcodes with more than
            10 adt UMIs is less than 100.")
  
  depth %>%
    filter(log_total_UMI_adt >= 0.2) -> sub_depth
    #filter(log_nfeat_gex >= 0.2) 
  dim(sub_depth)
  
  back <- filter(sub_depth, group == "background")
  flog.info("Barcode dimensions:", dim(back), capture = TRUE)
  
}

# ADT UMIs
back$log_total_UMI_adt[is.infinite(back$log_total_UMI_adt)] <- 0
max_adt <- mean(back$log_total_UMI_adt)+(2*sd(back$log_total_UMI_adt))

if(max_adt >= ccmin) {
  max_adt <- ccmin-0.2
}
min_adt <- mean(back$log_total_UMI_adt)-(3*sd(back$log_total_UMI_adt))

if(min_adt > 0.5) {
  min_adt <- 0.5
} else if(min_adt < 0) {
  min_adt <- 0.5
}

# GEX genes
back$log_nfeat_gex[is.infinite(back$log_nfeat_gex)] <- 0
head(back)
max_gen <- mean(back$log_nfeat_gex)+(2*sd(back$log_nfeat_gex))

if(max_gen > 80) {
  max_gen <- log10(80)
}

min_gen <- mean(back$log_nfeat_gex)-(3*sd(back$log_nfeat_gex))
flog.info("Mean(#genes) - 3sd(#genes):", min_gen, capture = TRUE)

if(min_gen > 0) {
  min_gen <- 0
} else if(min_gen < 0) {
  min_gen <- 0
}

gg_dsb_back <- ggplot(depth, aes(x = log_nfeat_gex, y = log_total_UMI_adt)) +
  geom_bin2d(bins = 100) +
  facet_wrap(~group) +
  ggtitle("Background vs cells based on ADT & GEX",
          subtitle = paste(table(depth$group), collapse = " ")) +
  geom_vline(xintercept = c(min_gen, max_gen), linetype = "dashed", 
             alpha = .4, color = "blue") +
  geom_hline(yintercept = c(min_adt, max_adt), linetype = "dashed", 
             alpha = .4, color = "blue") +
  xlim(0, NA) +
  ylim(0, NA) +
  scale_fill_gradient(low = "grey90", high = "darkred") +
  theme_bw() +
  xlab("# genes (GEX)") +
  ylab("# UMIs (ADT)")

# Filtering out lowly & highly seq background-barcodes -------------------------
back %>%
  filter(log_total_UMI_adt > min_adt &
           log_total_UMI_adt < max_adt &
           log_nfeat_gex > min_gen &
           log_nfeat_gex < max_gen) -> sub_back
  
flog.info("Final background dimensions:", dim(sub_back), capture = TRUE)

# Filtering out lowly & highly seq cell-barcodes -------------------------------
sub_cell <- filter(sub_depth, group == "cell") %>%
  filter(log_total_UMI_adt > max_adt) %>%
  filter(log_nfeat_gex >= log10(200)) %>%
  filter(log_nfeat_gex <= log10(cfmax))

# Assembling combined filtered background & cell seq metrics -------------------
sub_depth <- rbind(sub_back, sub_cell)

gg_dsb_sub_back <- ggplot(sub_depth, aes(x = log_nfeat_gex, 
                                     y = log_total_UMI_adt)) +
  geom_bin2d(bins = 100) +
  facet_wrap(~group) +
  scale_fill_gradient(low = "grey90", high = "darkred") +
  ggtitle("Background vs cells based on ADT & GEX (trimmed)",
          subtitle = paste(table(sub_depth$group), collapse = " ")) +
  xlim(0, NA) +
  ylim(0,NA) +
  theme_bw()

gg_sub_adt <- get_dens(sub_depth, v = "log_total_UMI_adt") +
  ggtitle("DSB inputs:",
          subtitle = paste(table(sub_depth$group),
                collapse = " ")) +
  scale_fill_manual(values = c("grey", "darkblue"))

write.table(sub_depth, gsub("\\.pdf", ".tsv", opt$outfile),
            sep = "\t", row.names = F, quote = F)

flog.info("Final barcodes:", table(sub_depth$group), capture = TRUE)

pdf(opt$outfile, width = 8, height = 14)

  plot_grid(
    plot_grid(plotlist = list(gg_adt_den,
                              gg_adt_dens,
                              gg_gex_den,
                              gg_gex_dens), 
              ncol = 2),
    gg_dsb_back, ncol = 1)
  
  plot_grid(plotlist = list(gg_dsb_sub_back,
                            gg_sub_adt),
            ncol = 1)
  
dev.off()

# If user provide limits for background and cell barcodes ----------------------

if(opt$bcmin != "None") {
  
  flog.info("Ploting background and cell containing barcodes based on \
            user's thresholds")
  
  gex_cell_barcodes <- filter(gex_depth,
                              nfeat > opt$cfmin &
                                nfeat < opt$cfmax)[["barcode_id"]]
  
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
