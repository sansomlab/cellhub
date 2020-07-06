library(tidyverse)
library(magrittr)
library(tidyverse)
library(colormap)
library(circlize)
library(ComplexHeatmap)
library(R.utils)
library(optparse)
options(stringsAsFactors = F)


option_list <- list(
  make_option(c("--basedir"), default="/well/singlecell/P190381/_zzz_/10X-count",
              help="the folder where project data of the demultiplexing. the name has to finish with _10x-count-_ to match with Santiago's new structure"),
  make_option(c("--samplename"), default="726070_GX06,726070_GX18",
              help="the name of the sample folders you want to skip to parse the project results (i.e. failed channels). Default to NULL"),
  make_option(c("--subset"), default="vireo-known_AF5e4,demuxlet-0.0001",
              help="the name of the subset of demultiplexing outputs you want specifically to plot, useful if you already know what options results you'd like to compare"),
  make_option(c("--outdir"), default="project.results",
              help="the name of the output folder")
)

message("Pipeline to parse data from genetic multiplexing outputs across channels, at the project level")  
opt <- parse_args(OptionParser(option_list=option_list))
if(is.null(opt$outdir)) { opt$outdir <- paste0(getwd(),"/")}

if(!grepl("\\/$", opt$outdir)){opt$outdir <- paste(opt$outdir, "/", sep = "")}
if(!grepl("\\/$", opt$basedir)){opt$basedir <- paste(opt$basedir, "/", sep = "")}

if (!file.exists(opt$outdir)){ dir.create(opt$outdir)}
print("Running with options:")
print(opt)


run<- opt$outdir
opt[which(opt=="NULL")] <- NULL
opt[which(opt=="None")] <- NULL



#skip samples instead of listing them

rdir<-system(paste0("ls -d ",opt$basedir,"/results.*"),intern=T)
gsub("results.","",basename(rdir)) ->rdir


if(!is.null(opt$samplename)){
  message("skipping these channels: ")
  opt$samplename<-strsplit(opt$samplename,",")[[1]]
  gsub(" ","",opt$samplename)->opt$samplename
  print(opt$samplename)
  rdir <- rdir[!rdir  %in% opt$samplename]
}

if(!is.null(opt$subset)){
  opt$subset<-strsplit(opt$subset,",")[[1]]
  gsub(" ","",opt$subset)->opt$subset
}else{
  cnames<-read.table( paste0(opt$basedir, "results.",rdir[1],"/",rdir[1],"_SingleCellMetadata_demultiplexing_results.tsv.gz"), header=T, sep="\t", check.names = F)
  opt$subset <- colnames(cnames)[colnames(cnames)!="BARCODE"]
}



dlist<- list()
for (sub in opt$subset){
  dlist[[sub]]<- list()
  for (sam in rdir){
    dlist[[sub]][[sam]] <- read.table( paste0(opt$basedir, "results.",sam,"/",sam,"_SingleCellMetadata_demultiplexing_results.tsv.gz"), header=T, sep="\t", check.names = F)[,c("BARCODE",sub)]
    #dlist[[sub]][[sam]][,"BARCODE"]<-gsub("-1",paste0("-",sam),dlist[[sub]][[sam]][,"BARCODE"])
    dlist[[sub]][[sam]][,"sample"]<-rep(sam)
  }
}

res.list<-lapply(dlist, function(x) do.call("rbind",x))
test<-lapply(res.list, function(x) mutate(x, barcode_id=paste(BARCODE, sample, sep="-")) %>% select(-c(BARCODE,sample) )) 
data.write <-Reduce(function(x,y) merge(x,y,by="barcode_id"), test) 
write.table(data.write, file=paste0(run, "all_demultiplexing_results.tsv"), row.names = F, col.names = T, sep="\t", quote = F)
if(!file.exists(paste0(run,"comparison_selection"))) {dir.create(paste0(run,"comparison_selection"))}

for(A in opt$subset) {
  for(B in opt$subset) {
    if(A!=B){
      unclass(table(x=data.write[,A], y=data.write[,B])) ->mat
    
      Heatmap(mat, cluster_rows = F,cluster_columns = F,
              column_title =B,
              row_title = A,
              row_names_side = "left", 
              column_title_side = "bottom",
              row_names_gp = gpar(fontsize=9),
              column_names_gp = gpar(fontsize=9),
              cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 8))
            }) ->g
    pdf(paste0(run,"comparison_selection/",A,"_VS_",B,".pdf"), width = 8, height = 8)
    ComplexHeatmap::draw(g)
    dev.off()
    }
  }
}

data.write ->df
as.character(unique(unlist(df[,opt$subset])))->gnr
c("DOUBLET","UNASSIGNED",gnr[!gnr %in% c("DOUBLET","UNASSIGNED")])->gnr
sum(!gnr %in% c("DOUBLET","UNASSIGNED")) ->vc
colormap(colormap=colormaps$portland,vc) ->ccol
c("#7d00b3ff","#93ff00ff", ccol) ->vcol

df %>%  
  pivot_longer( names_to="method", values_to="general", cols=opt$subset) %>% 
  mutate(general=factor(general, levels=gnr),
          sample=sapply(strsplit(barcode_id,"-"),'[[',3) ) %>% 
  group_by(sample,general, method) %>% 
  summarise(tot=n()) %>% 
  group_by(sample,method) %>% 
  mutate(persample=sum(tot), 
         percent=100*tot/persample) %>% 
  ggplot(aes(method, percent, fill=general)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values = vcol, 
                    breaks=gnr, 
                    limits=gnr) +
  theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1), 
        legend.position="top", legend.key.size = unit(0.3, "cm")) +
        guides(colour = guide_legend(nrow = 3))+
  facet_grid(~sample, scales="free_x") ->g

ggsave(g, filename = paste0(run,"methods.output.pdf"), width = 13, height = 6)




for (sub in opt$subset){
  write.table(res.list[[sub]],file=paste0(run, sub, "_demultiplexing_results.tsv"), row.names = F, col.names = T, sep="\t", quote = F)
  as.character(unique(unlist(res.list[[sub]][,sub])) )->gnr
  c("DOUBLET","UNASSIGNED",gnr[!gnr %in% c("DOUBLET","UNASSIGNED")])->gnr
  sum(!gnr %in% c("DOUBLET","UNASSIGNED")) ->vc
  colormap(colormap=colormaps$portland,vc) ->ccol
  c("#7d00b3ff","#93ff00ff", ccol) ->vcol
  
  res.list[[sub]] %>% 
    rename(general=sub) %>% 
    mutate(general=factor(general, levels=gnr)) %>% 
    group_by(sample,general) %>% 
    summarise(tot=n()) %>% 
    group_by(sample) %>% 
    mutate(persample=sum(tot), 
           percent=100*tot/persample) %>% 
    ggplot(aes(sample, percent, fill=general)) +
    geom_bar(stat="identity", position="stack", color="black") +
    scale_fill_manual(values = vcol, 
                      breaks=gnr, 
                      limits=gnr) +
    theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1),
          legend.position="top", legend.key.size = unit(0.3, "cm")) +
    guides(colour = guide_legend(nrow = 3))+
      ggtitle(sub)->g1
  
  ggsave(g1, filename = file.path(paste0(run, sub, "_demultiplexing_results.pdf")), width=10, height = 7)
  
  res.list[[sub]] %>% 
    rename(general=sub) %>% 
    mutate(general=factor(general, levels=gnr)) %>% 
    group_by(sample,general) %>% 
    summarise(tot=n()) %>% 
    group_by(sample) %>% 
    mutate(persample=sum(tot), 
           percent=100*tot/persample) %>% 
    ggplot(aes(sample, tot, fill=general)) +
    geom_bar(stat="identity", position="stack", color="black") +
    scale_fill_manual(values = vcol, 
                      breaks=gnr, 
                      limits=gnr) +
    theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1),
          legend.position="top", legend.key.size = unit(0.3, "cm")) +
    guides(colour = guide_legend(nrow = 3))+
    ggtitle(sub) + ylab("cell yield") ->g1
  
  ggsave(g1, filename = file.path(paste0(run, sub, "_demultiplexing_results_cellcount.pdf")), width=10, height = 7)
  
  
  
}


message("donut")