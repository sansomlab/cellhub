library(ggplot2)

pca_var = read.table(gzfile("pca_variance_ratios.tsv.gz"), header=T)

pca_var$nPC <- as.numeric(rownames(pca_var))
colnames(pca_var) <- c("var", "nPC")
gp <- ggplot(data=pca_var, aes(x=nPC, y=var)) + geom_point(size=0.7)
gp <- gp + theme_bw()
gp <- gp + scale_y_log10() + ylab("variance ratio")
gp <- gp + scale_x_continuous(breaks = seq(0, 100, 5)) + theme(panel.grid.minor =   element_blank())

ggsave(filename = "figures.dir/pca_variance_ratio_ggplot.pdf", plot = gp, width = 6, height = 4)
#options(bitmapType = 'cairo', device = 'png')
#ggsave(filename = "figures.dir/pca_variance_ratio_ggplot.png", plot = gp, width = 6, type="cairo", height = 4)