#
# Quality Control of a 10X GEX mapping output
#

" Quality Control of a 10X GEX mapping output

Usage: build_qc_mapping_reports.R --tenxfolder=<path> --sample_id=<string> --specie=<value> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01.
  --tenxfolder=<path>   10X GEX output folder, /outs parent.
  --sample_id=<string> Sample/channel id to call output report
  --specie=<value>     Either 'hg' or 'mm'.
  --outfolder=<file>   Path to results folder.

Author: 
  Cesar Prada
e-mail: 
  cesarprada[at]email.com

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'batchelor', 'ggplot2',
          'cowplot', "dplyr", "ggwordcloud")

if (!require("BiocManager", character.only = TRUE)) {
	install.packages("BiocManager")
	BiocManager::install()
} else {
	ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
	if (any(!ipkgs)) {
		BiocManager::install(pkgs[!ipkgs])
	} else {
		message("\n\nCool! your machine has everything is needed.\n\n")
	}
}

# -------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(ggwordcloud))
# ------------------------------

# --------------
# Main function
# --------------

plot_qc <- function(cds, mt_genes, rb_genes, ds, jobname) {

    # ---------------------------
    # Feature centered histograms
    # ---------------------------

    logcount <- log10(rowSums(cds[["exprs"]]))
    features <- data.frame('feat_count_depth' = rowSums(cds[["exprs"]]),
                           'feat_log10_count_depth' = ifelse(is.infinite(logcount),
						       0.001, logcount),
			   'feat_cell_depth' = rowSums(cds[["exprs"]] != 0),
			   'gene_short_name' = cds[["fData"]][rownames(cds[["exprs"]]), ]$gene_short_name
			   )

    features[['id']] <- rownames(features)

    librarysize <- sum(features$feat_count_depth)
    features$pctlib <- features$feat_count_depth/librarysize
    
    libsizecell <- colSums(cds[["exprs"]])

    top20feat <- head(arrange(features, -feat_count_depth),20)$id
    top20exp <- cds[["exprs"]][top20feat, ]
    pctlibsizecell <- lapply(seq(ncol(top20exp)), function(i) top20exp[, i]/libsizecell[i])
    names(pctlibsizecell) <- colnames(top20exp)
    
    top20 <- do.call(rbind, 
	      lapply(names(pctlibsizecell), function(cell) {
		  data.frame('pctlibcellsize'=pctlibsizecell[[cell]],
			     'id'=names(pctlibsizecell[[cell]]),
			     'cell'=cell)
		  })
	     )

    print(head(top20))
    top20 <- merge(top20, select(features, c(`id`, gene_short_name)), by='id')
    print("Done")
    
    gtop20 <- ggplot(top20, aes(x = gene_short_name, y=pctlibcellsize)) +
	    geom_jitter(aes(group = gene_short_name), alpha=.2, color = "grey", shape=".") +
	    geom_violin(alpha=.3, fill = "grey") +
	    scale_x_discrete(limits=rev(head(arrange(features, -feat_count_depth),20)$gene_short_name)) +
	    coord_flip() +
	    theme_classic() +
	    theme(axis.text=element_text(size=5)) +
	    ggtitle(paste(ds, "Top 20 deepest genes"),
		    subtitle = paste0("They account for ~",
				     ceiling(sum(features[top20feat,]$pctlib)*100),
				     "% of total-counts"))

     g1 <- ggplot(`features`, aes(`feat_count_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_count_depth)/100)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_count_depth))), 
		 x = mean(features$feat_count_depth), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_count_depth, .99)), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature count depth"))

    g1_1 <- ggplot(features, aes(feat_count_depth)) +
	  geom_histogram(color = "grey", fill = "grey", 
			 bins = (ceiling(max(features$feat_count_depth)/100))) +
	  theme_classic() +
	  xlim(c(-0.5, ceiling(quantile(features$feat_count_depth, 0.5)))) +
	  ggtitle(paste(ds, "Feature count depth"),
		  subtitle = "Quantile 50")

    g2_wc <- ggplot(head(arrange(features, -feat_count_depth), nrow(features)*0.01), 
		    aes(label = gene_short_name, 
				size = feat_count_depth, 
				color = feat_count_depth)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% count-deepest features")


    g2 <- ggplot(features, aes(feat_log10_count_depth)) +
	geom_histogram(color = "grey", fill = "grey", bins=100) +
	geom_vline(xintercept = round(quantile(features$feat_log10_count_depth, 0.5), digits = 2),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       round(mean(features$feat_log10_count_depth), digits=2)), 
		 x = mean(features$feat_log10_count_depth), 
		 y = max(hist(features$feat_log10_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_log10_count_depth),
					   by = max(features$feat_log10_count_depth)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "log10(feature count depth)"))

    g1 <- ggplot(`features`, aes(`feat_count_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_count_depth)/200)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_count_depth))), 
		 x = mean(features$feat_count_depth), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/200)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_count_depth, .99)), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature count depth"))


    g2_1 <- ggplot(features, aes(feat_log10_count_depth)) +
	  geom_histogram(color = "grey", fill = "grey", binwidth = 0.01) +
          xlim(c(-.1, quantile(features$feat_log10_count_depth, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "log10(feature count depth)"),
		  subtitle = "Quantile 50")

    g2_2 <- ggplot(`features`, aes(`feat_cell_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_cell_depth)/100)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_cell_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_cell_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_cell_depth))), 
		 x = mean(features$feat_cell_depth), 
		 y = max(hist(features$feat_cell_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_cell_depth),
					   by = max(features$feat_cell_depth)/ceiling(max(features$feat_cell_depth)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_cell_depth, .99)), 
		 y = max(hist(features$feat_cell_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_cell_depth),
					   by = max(features$feat_cell_depth)/ceiling(max(features$feat_cell_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature cell depth"))

    g2_2_wc <- ggplot(head(arrange(features, -feat_cell_depth), nrow(features)*0.01), 
		    aes(label = gene_short_name, 
				size = feat_cell_depth, 
				color = feat_cell_depth)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% cell-deepest features")



    # ------------------
    # Cell based density
    # ------------------

    logcount <- log10(colSums(cds[["exprs"]]))
    cell <- data.frame('cell_depth_counts' = colSums(cds[["exprs"]]),
                       'cell_log10_depth_counts' = ifelse(is.infinite(logcount),
							  0.1, logcount),
		       'cell_depth_genes_detected' = colSums(cds[["exprs"]] != 0)
		       )

    cell[['barcode']] <- rownames(cell)
    print(head(cell))

    # ---
    # Mitochondrial gene information
    # ---

    if(length(mt_genes) == 0) {
	    mtc <- setNames(rep(0, ncol(cds[["exprs"]])), colnames(cds[["exprs"]]))
    } else {
	    mtc <- colSums(cds[["exprs"]][mt_genes, ])
    }
    mt <- data.frame(mt_counts = mtc)
    mt$pct_mt_counts <- mt$mt_counts/cell$cell_depth_counts
    mt$barcode <- rownames(mt)
    print(head(mt))
    cell <- merge(cell, mt, by = "barcode")

    print(head(cell))

    # ---
    # Ribosomal
    # ---

    if(length(rb_genes) == 0) {
	    rbc <- setNames(rep(0, ncol(cds[["exprs"]])), colnames(cds[["exprs"]]))
    } else {
	    rbc <- colSums(cds[["exprs"]][rb_genes, ])
    }
    rb <- data.frame(rb_counts = rbc)
    rb$pct_rb_counts <- rb$rb_counts/cell$cell_depth_counts
    rb$barcode <- rownames(rb)

    print("Rb genes")
    print(head(cell))
    cell <- merge(cell, rb, by = "barcode")
    print(head(cell))

    cell %>%
	dplyr::arrange(-cell_depth_counts) %>%
	mutate(`rank` = seq(nrow(cell))) -> cell
    print(head(cell,2))
   
    rownames(cell) <- cell$barcode

    g3 <- ggplot(cell, aes(cell_depth_counts)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 1000) +
	geom_vline(xintercept = ceiling(quantile(cell$cell_depth_counts, 0.5)),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(cell$cell_depth_counts))), 
		 x = mean(cell$cell_depth_counts), 
		 y = max(hist(cell$cell_depth_counts, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_depth_counts),
					   by = max(cell$cell_depth_counts)/1000),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Count depth"))

    g3_1 <- ggplot(cell, aes(cell_depth_counts)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
	  xlim(c(-.5, ceiling(quantile(cell$cell_depth_counts, 0.5)))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Count depth - Quantile 50"))

    g4 <- ggplot(cell, aes(cell_log10_depth_counts)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell$cell_log10_depth_counts, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       round(mean(cell$cell_log10_depth_counts), digits=2)), 
		 x = mean(cell$cell_log10_depth_counts), 
		 y = max(hist(cell$cell_log10_depth_counts, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_log10_depth_counts),
					   by = max(cell$cell_log10_depth_counts)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "log10(count depth)"))

    g4_1 <- ggplot(cell, aes(cell_log10_depth_counts)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 50) +
	  xlim(c(-.1, quantile(cell$cell_log10_depth_counts, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "log10(count depth) - Quantile 50"))

    g5 <- ggplot(cell, aes(cell_depth_genes_detected)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell$cell_depth_genes_detected, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(cell$cell_depth_genes_detected))), 
		 x = mean(cell$cell_depth_genes_detected), 
		 y = max(hist(cell$cell_depth_genes_detected, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_depth_genes_detected),
					   by = max(cell$cell_depth_genes_detected)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Gene depth"))

    g5_1 <- ggplot(cell, aes(cell_depth_genes_detected)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", binwidth = 1) +
	  xlim(c(-.1, quantile(cell$cell_depth_genes_detected, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Gene depth - Quantile 50"))

    g6 <- ggplot(cell, aes(x = `rank`, y = cell_depth_counts)) +
	geom_point(aes(color = pct_mt_counts), 
		   stat = "identity", alpha = .2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	ylim(c(0, max(cell$cell_depth_counts))) +
	theme_classic() +
	ggtitle(paste(ds, "Depth vs. rank"))

    g6_1 <- ggplot(cell, aes(x = `rank`, y = cell_depth_counts)) +
	  geom_point(aes(color = pct_rb_counts), 
		     stat = "identity", alpha = .2, size = 0.5) +
	  scale_color_gradient(low = "pink",  high = "black") +
	  ylim(c(0, max(cell$cell_depth_counts))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Depth vs. rank"))

    g7 <- ggplot(cell, aes(x = cell_depth_counts, y = cell_depth_genes_detected)) +
	geom_point(aes(color = pct_mt_counts), stat = "identity", 
		   alpha = 0.2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	theme_classic() +
	ggtitle(paste(ds, "Genes-detected vs. depth"),
		subtitle = "% Mitochondrial-counts")

    g8 <- ggplot(cell, aes(x = cell_depth_counts, y = cell_depth_genes_detected)) +
	geom_point(aes(color = pct_rb_counts, fill = pct_rb_counts), 
		   stat = "identity", alpha = 0.2, size = 0.5) +
	scale_color_gradient(high = "black",  low = "pink") +
	scale_fill_gradient(high = "black",  low = "pink") +
	theme_classic() +
	ggtitle(paste(ds, "Genes-detected vs. depth"),
	       subtitle = "% Ribosomal-counts")

    if(max(cell$pct_mt_counts) == 0) {


	    g9 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_mt_counts)) +
		geom_point(color = "red", alpha = .2, size = 0.5) +
		theme_classic() +
		ggtitle(paste(ds, "% Mitochondrial-counts vs. depth"))
	

	    g10 <-  ggplot(cell, aes(pct_mt_counts)) +
		        geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_mt_counts),
				   linetype = "dashed", color = "black") + 
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Mitochondrial-counts density"))

    } else {

	    g9 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_mt_counts)) +
		    geom_point(color = "red", alpha = .2, size = 0.5) +
		    geom_density_2d(color = "black", alpha = 0.5) +
		    theme_classic() +
		    ggtitle(paste(ds, "% Mitochondrial-counts vs. depth"))

	    g10 <- ggplot(cell, aes(pct_mt_counts)) +
			geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_mt_counts),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("mean ~", 
							round(mean(cell$pct_mt_counts), digits=2)), 
					  	 x = mean(cell$pct_mt_counts), 
						 y = max(hist(cell$pct_mt_counts, 
							      breaks = seq(from = 0, 
									   to = max(cell$pct_mt_counts),
									   by = max(cell$pct_mt_counts)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Mitochondrial-counts density"))
    }

    if(max(cell$pct_rb_counts) == 0) {

	    g11 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_rb_counts)) +
		geom_point(alpha = 0.4, color = "pink", size = 0.5) +
		theme_classic() +
		ggtitle(paste(ds, "% Ribosomal-counts vs. depth"))


	    g12 <-  ggplot(cell, aes(pct_rb_counts)) +
		        geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_rb_counts),
				   linetype = "dashed", color = "black") + 
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Ribosomal-counts density"))

    } else {

	    g11 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_rb_counts)) +
		geom_point(alpha = 0.4, color = "pink", size = 0.5) +
		geom_density_2d(color = "black", alpha = 0.6) +
		theme_classic() +
		ggtitle(paste(ds, "% Ribosomal-counts vs. depth"))

	    g12 <- ggplot(cell, aes(pct_rb_counts)) +
			geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_rb_counts),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("mean ~", 
							round(mean(cell$pct_rb_counts), digits=2)), 
					  	 x = mean(cell$pct_rb_counts), 
						 y = max(hist(cell$pct_rb_counts, 
							      breaks = seq(from = 0, 
									   to = max(cell$pct_rb_counts),
									   by = max(cell$pct_rb_counts)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Ribosomal-counts density"))

    }

    return(list("ps1" = list(g1, g2_wc, g2, gtop20, g2_2, g2_2_wc, g3, g3_1, g4, g4_1),
		"ps2" = list(g5, g5_1, g6, g7, g9, g10, g6_1, g8, g11, g12)))
}

# Twiking monocle3 load_cellranger_data function ()
# to return a list rather than a cds object

load_cellranger_data <- function(pipestance_path=NULL, genome=NULL,
                                 barcode_filtered=TRUE, umi_cutoff = 100) {
  # check for correct directory structure
  if (!dir.exists(pipestance_path))
    stop("Could not find the pipestance path: '", pipestance_path,"'.
         Please double-check if the directory exists.\n")

  od = file.path(pipestance_path)

  if (!dir.exists(od))
    stop("Could not find the pipestance output directory: '",
         file.path(pipestance_path),
         "'. Please double-check if the directory exists.\n")

  v3p = file.path(od)
  v2p = file.path(od, "filtered_gene_bc_matrices")
  v3d = dir.exists(v3p)
  
  matrix_dir = ifelse(v3d, v3p, v2p)
  print(matrix_dir)

  if(!dir.exists(matrix_dir))
    stop("Could not find directory: ", matrix_dir)

  if(v3d) {
    features.loc <- file.path(matrix_dir, "features.tsv.gz")
    barcode.loc <- file.path(matrix_dir, "barcodes.tsv.gz")
    matrix.loc <- file.path(matrix_dir, "matrix.mtx.gz")
  } else {
    barcode.loc <- file.path(matrix_dir, genome, "barcodes.tsv")
    features.loc <- file.path(matrix_dir, genome, "genes.tsv")
    matrix.loc <- file.path(matrix_dir, genome, "matrix.mtx")
  }
  if (!file.exists(barcode.loc)){
    stop("Barcode file missing")
  }
  if (!file.exists(features.loc)){
    stop("Gene name or features file missing")
  }
  if (!file.exists(matrix.loc)){
    stop("Expression matrix file missing")
  }
  data <- Matrix::readMM(matrix.loc)

  feature.names = utils::read.delim(features.loc,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  # Duplicate row names not allowed
  feature.names$V1 = make.unique(feature.names$V1)
  if(dim(data)[1] != length(feature.names[,1])) {
    stop(sprintf(paste("Mismatch dimension between gene file: \n\t %s\n and",
                       "matrix file: \n\t %s\n"), features.loc, matrix.loc))
  }
  if(v3d) {
    # We will only load GEX data for the relevant genome
    data_types = factor(feature.names$V3)
    print(data_types)
    allowed = data_types %in% c("Gene Expression", "Antibody Capture")
    if(!is.null(genome)) {
      # If not multigenome, no prefix will be added and we won't filter out
      # the one genome
      gfilter = grepl(genome, feature.names$V1)
      if(any(gfilter)) {
        allowed = allowed & grepl(genome, feature.names$V1)
      } else {
        message(paste("Data does not appear to be from a multi-genome sample,",
                      "simply returning all gene feature data without",
                      "filtering by genome."))
      }

    }
    data = data[allowed, ]
    feature.names = feature.names[allowed,1:2]
  }
  colnames(feature.names) = c("id", "gene_short_name")
  rownames(data) = feature.names[,"id"]
  rownames(feature.names) = feature.names[,"id"]

  barcodes <- utils::read.delim(barcode.loc, stringsAsFactors=FALSE, header=FALSE)
  if (dim(data)[2] != length(barcodes[,1])) {
    stop(sprintf(paste("Mismatch dimension between barcode file: \n\t %s\n",
                       "and matrix file: \n\t %s\n"), barcode.loc,matrix.loc))
  }
  barcodes$V1 = make.unique(barcodes$V1)
  colnames(data) = barcodes[,1]
  pd = data.frame(barcode=barcodes[,1], row.names=barcodes[,1])
  data <- data[,Matrix::colSums(data) > umi_cutoff]
  pd <- pd[colnames(data),, drop=FALSE]

  gbm <- list("exprs" = data,
	      "pData" = pd,
	      "fData" = feature.names)

  return(gbm)
}


print("Functions loaded")

# ---------------
# --- Parameters
# ---------------

inputfolder <- arguments$tenxfolder
sp <- arguments$specie
outputfolder <- arguments$outfolder
jobname <- arguments$sample_id

print("Parameters loaded")

# Read single_cell_object
cds <- load_cellranger_data(inputfolder, umi_cutoff = 0)

if(class(cds) != "list") {

    stop(paste("Sorry", print(class(cds)), "objects are not supported!",
	       "Try cell_data_set (Monocle3) instead!"))

}

print(str(cds))

if(sp == "mm") {
	
	if(colnames(cds[["fData"]])[1] == "gene_id") {

		colnames(cds[["fData"]])[1] <- "id"
	}

	# Mitochondrial genes
	mt_genes <- cds[["fData"]][grepl("^mt-", cds[["fData"]]$gene_short_name), ]$id
	message("Mitochondrial genes detected:")
	print(mt_genes)
	if(length(mt_genes) == 0) {
		message("No mitochondrial genes detected") 
	}

	# Ribosomal genes
	rb_genes <- cds[["fData"]][grepl("^Rps|^Rpl", cds[["fData"]]$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(rb_genes)
	if(length(rb_genes) == 0) {
		message("Ribosomal genes detected")
	}

} else if (sp == "hg") {
	
	# Mitochondrial genes
	mt_genes <- cds[["fData"]][grepl("^MT-", cds[["fData"]]$gene_short_name), ]$id
	message("Mitochondrial genes detected:")
	print(mt_genes)

	if(length(mt_genes) == 0) {
		message("No mitochondrial genes detected") 
	}
	
	# Ribosomal genes
	rb_genes <- cds[["fData"]][grepl("^RPS|^RPL", cds[["fData"]]$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(rb_genes)
	if(length(rb_genes) == 0) {
		message("Ribosomal genes detected")
	}

} else {

	stop("Please provide a valid specie parameter")

}

# --- Run
res <- plot_qc(cds, mt_genes, rb_genes, jobname)

# --- Write QC plots
title <- ggdraw() + 
  draw_label(
    jobname,
    fontface = 'bold',
    #x = 0,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )

qc_file <- paste0(outputfolder, '/', jobname, '_', 'qcmetrics_report.pdf')

pdf(qc_file, height = 11.69, width = 8.27)
    plot_grid(title, plot_grid(plotlist = res$ps1, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
    plot_grid(title, plot_grid(plotlist = res$ps2, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

message(paste("The QC report ", qc_file, "has been created."))
