library(data.table)
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
  make_option(c("--demultiplexing"), default="demuxlet,vireo",
              help="the folder where project data of the demultiplexing"),
  make_option(c("--samplename"), default="726070_GX06",
              help="the name of the sample folder"),
  make_option(c("--subset"), default="vireo.726070_GX06-known_AF5e4,demuxlet.726070_GX06-0.0001",
              help="the name of the subset of demultiplexing outputs you want specifically to plot, useful if you already know what options results you'd like to compare"),
  make_option(c("--outdir"), default=NULL,
              help="the name of the sample folder")
)
  
message("Pipeline to parse data from genetic multiplexing outputs")  
opt <- parse_args(OptionParser(option_list=option_list))
if(is.null(opt$outdir)) { opt$outdir <- paste0(getwd(),"/")}
if(!grepl("\\/$", opt$outdir)){opt$outdir <- paste(opt$outdir, "/", sep = "")}

print("Running with options:")
print(opt)

#parse adn format outdir, print opts statement
run<- opt$outdir
opt[which(opt=="NULL")] <- NULL


opt$demultiplexing<-strsplit(opt$demultiplexing,",")[[1]]
if(length(opt$demultiplexing)<2){ #uniform names results
  if (opt$demultiplexing=="demuxlet"){ 
    dend<- ".best"
  }else if(opt$demultiplexing=="vireo") {
    dend<- "donor_ids.tsv"
  }else{
    stop("can't recognise this multiplexing output")
  }
  pstring <- paste0( opt$basedir,"-", opt$demultiplexing, "/",opt$samplename,"*/*",dend)
  ll <- system(paste0("ls ",pstring), intern=T)
  dlist <-list()
  if(opt$demultiplexing=="demuxlet") {
    cid<-c("BARCODE","BEST")
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("demuxlet.",basename(dirname(ll[idx])) )
      general<-sapply(strsplit(dlist[[idx]][,2],"-"),"[[",1)
      specific<-sapply(strsplit(dlist[[idx]][,2],"-"),"[[",2)
      which(general=="SNG") ->nex
      general[nex]<-specific[nex]
      which(general == "AMB") ->NEX
      general[NEX] <- "UNASSIGNED"
      which(general == "DBL") ->DEX
      general[DEX] <- "DOUBLET"
      dlist[[idx]][,2] <- general
      nex <- NEX <-DEX <- NULL
    }
    names(dlist) <- basename(dirname(ll)) 
  }
  if(opt$demultiplexing=="vireo"){
    cid<-c("cell","donor_id")
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("vireo.",basename(dirname(ll[idx]))) 
      colnames(dlist[[idx]])[1] <-"BARCODE"
      dlist[[idx]][,2]<-toupper(dlist[[idx]][,2])
    }

    names(dlist) <- basename(dirname(ll)) 
  }
  data.use<- Reduce(function(x,y) merge(x,y,by="BARCODE"), dlist)
  message("-----------------")
  message("cells demultiplexed :")
  print(nrow(data.use))
  
}else{
  message("multiple demultiplexing detected, parsing all")
  dm<-opt$demultiplexing[which(opt$demultiplexing=="demuxlet")]
  dstring<- paste0( opt$basedir,"-demuxlet/",opt$samplename,"*/*.best")
  vm<-opt$demultiplexing[which(opt$demultiplexing=="vireo")]
  vstring<- paste0( opt$basedir, "-vireo/",opt$samplename,"*/donor_ids.tsv")
  
  demtemp <- "demuxlet"
  cid<-c("BARCODE","BEST")
  ll <- system(paste0("ls ",dstring), intern=T)
  dlist <-list()
  for( idx in 1:length(ll)){
    if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
    dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
    colnames(dlist[[idx]])[2] <- paste0("demuxlet.",basename(dirname(ll[idx])) )
    general<-sapply(strsplit(dlist[[idx]][,2],"-"),"[[",1)
    specific<-sapply(strsplit(dlist[[idx]][,2],"-"),"[[",2)
    which(general=="SNG") ->nex
    general[nex]<-specific[nex]
    which(general == "AMB") ->NEX
    general[NEX] <- "UNASSIGNED"
    which(general == "DBL") ->DEX
    general[DEX] <- "DOUBLET"
    dlist[[idx]][,2] <- general
    nex <- NEX <-DEX <- NULL
  }
  names(dlist) <- basename(dirname(ll)) 
  data.demuxlet<- Reduce(function(x,y) merge(x,y,by="BARCODE"), dlist)
  
  demtemp <- "vireo"
  cid<-c("cell","donor_id")
  ll <- system(paste0("ls ",vstring), intern=T)
  dlist <-list()
  for( idx in 1:length(ll)){
    if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
    dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
    colnames(dlist[[idx]])[2] <- paste0("vireo.",basename(dirname(ll[idx]))) 
    colnames(dlist[[idx]])[1] <-"BARCODE"
    dlist[[idx]][,2]<-toupper(dlist[[idx]][,2])
  }
  names(dlist) <- basename(dirname(ll)) 
  data.vireo<- Reduce(function(x,y) merge(x,y,by="BARCODE"), dlist)
  message("-----------------")
  message("cells used for vireo :")
  print(nrow(data.vireo))
  message(" ")
  message("cells used for demuxlet :")
  print(nrow(data.demuxlet))
  message("-----------------")
  data.use <- merge(data.vireo, data.demuxlet, by="BARCODE")
  message("Merged cells demultiplexed :")
  print(nrow(data.use))
}


data.use %>% reshape2::melt(id.vars="BARCODE") %>% 
  group_by(variable, value) %>% summarise(tot=n()) %>% 
  group_by(variable) %>% mutate(persample=sum(tot), percent=100*tot/persample) -> infomat

infomat %>% 
  pivot_wider(id_cols=value,
              values_from = percent, 
              names_from = variable) %>% 
  replace(is.na(.),0) ->save.sng

infomat %>% mutate(General=ifelse(!(value %in% c("UNASSIGNED","DOUBLET")), "SINGLET",value )) %>% 
  group_by(variable, General) %>% summarise(tot=sum(tot)) %>%
  group_by(variable) %>% mutate(persample=sum(tot), percent=100*tot/persample) -> infomat.general

infomat.general %>% 
  pivot_wider(id_cols=General,
              values_from = percent, 
              names_from = variable) %>% 
  replace(is.na(.),0) ->save.general

output_path <-  file.path(paste0(run, opt$samplename, "_SingleCellMetadata_demultiplexing_results.txt"))
write.table(data.use,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)

output_path <-  file.path(paste0(run, opt$samplename, "_General_demultiplexing_results.txt"))
write.table(save.general,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)

output_path <-  file.path(paste0(run, opt$samplename, "_Singlets_demultiplexing_results.txt"))
write.table(save.sng,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)

c("#4e79e9ff","#93ff00ff","#7d00b3ff") ->scol


infomat.general %>% 
  mutate(General=factor(General, levels=c("UNASSIGNED","DOUBLET","SINGLET"))) %>% 
  ggplot(aes(variable, percent, fill=General)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values = scol, 
                    breaks=c("DOUBLET","UNASSIGNED","SINGLET"), 
                    limits=c("DOUBLET","UNASSIGNED","SINGLET")) +
  theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1),
        legend.position="top", legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))->g

ggsave(g, filename = file.path(paste0(run, opt$samplename, "_General_demultiplexing_results.pdf")), width=10, height = 7)

unique(unlist(infomat$value)) ->gnr
gnr <- c("UNASSIGNED",gnr[!gnr=="UNASSIGNED"])
sum(!gnr %in% c("DOUBLET","UNASSIGNED")) ->vc
colormap(colormap=colormaps$portland,vc) ->ccol
c("#93ff00ff","#4e79e9ff", ccol) ->vcol


infomat %>% mutate(value=factor(value, levels = gnr)) %>% 
  ggplot(aes(variable, percent, fill=value)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values = vcol, 
                    breaks=gnr, 
                    limits=gnr) +
  theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1),
        legend.position="top", legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3))->g1

ggsave(g1, filename = file.path(paste0(run, opt$samplename, "_Singlets_demultiplexing_results.pdf")), width=10, height = 7)

colnames(data.use)[!colnames(data.use)=="BARCODE"]->cols

ddf<-matrix(0,nrow=length(cols), ncol = length(cols))
rownames(ddf)<- colnames(ddf)<- cols
col_fun = colorRamp2(c(0, 50, 100), c("dodgerblue4", 'peachpuff', 'deeppink4'))

for(a in cols){
  for(b in cols){
    ddf[a,b] <- mclust::adjustedRandIndex(data.use[,a],data.use[,b])  
  }
}
  
Heatmap(ddf, cluster_rows = F,cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.4f", ddf[i, j]), x, y, gp = gpar(fontsize = 8))
        }) ->gg

pdf(paste0(run, opt$samplename,"_adjustedRand_scores.pdf"), 10,10)
gg
dev.off()
#draw gg should be the function but for some reason it doesn't take it! 

ddf %>% as.data.frame() %>% rownames_to_column("Var1") %>% 
  pivot_longer(cols, names_to = "Var2") %>% 
  filter(Var2 >Var1) %>% 
  arrange(-value) %>% rename(AdjustedRandIndex=value)->arranged.rand

output_path <-  file.path(paste0(run, opt$samplename, "_RandIndexScores.txt"))
write.table(arranged.rand,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)

unique(infomat.general$persample) ->fl

infomat.general %>% 
  ungroup() %>% 
  mutate(General=factor(General)) %>% 
  complete(General,variable, 
           fill=list(tot=0, persample=fl, percent=0)) %>% 
  filter(General=="SINGLET") %>% 
  top_n(percent,n=5) ->best



infomat.general %>% 
  ungroup() %>% 
  mutate(General=factor(General)) %>% 
  complete(General,variable, 
           fill=list(tot=0, persample=fl, percent=0)) %>% 
filter(General=="SINGLET") %>% 
  top_n(percent,n=-5) ->worst  


output_path <-  file.path(paste0(run, opt$samplename, "_top5_demultiplexing_results.txt"))
write.table(best,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)

output_path <-  file.path(paste0(run, opt$samplename, "_worst5_demultiplexing_results.txt"))
write.table(worst,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)


if (!is.null(opt$subset)) {
  opt$subset <-strsplit(opt$subset,",")[[1]]
  data.use %>% select(ends_with(c(opt$subset, as.character(best$variable) ))) ->data.subset
}else {
  data.use %>% select(ends_with(as.character(best$variable) )) ->data.subset
}


dir.create(paste0(run,"comparison_selection"))
cols <- colnames(data.subset)
for(a in cols){
  for(b in cols){
    n<-colSums(table(data.subset[,a],data.subset[,b])) #colsums over B
    norm.mat<-t(apply(table(data.subset[,a],data.subset[,b]),1, function(x) 100*x/n))
    Heatmap(norm.mat, 
            column_title = b,
            row_title = a,
            column_title_side = "bottom",
            row_title_side = "left",
            #cluster_rows = F,
            col = col_fun, 
            show_column_dend = T, 
            column_names_side = "top",
            row_names_side = "left", 
            rect_gp = gpar(col = "grey8", lwd = 0.6),
            name = paste0("Percent\nover\n",b)) ->gmat
    aa <- gsub(paste0(".",opt$samplename),"",a)
    bb <- gsub(paste0(".",opt$samplename),"",b)
    pdf(paste0(run,"comparison_selection/",opt$samplename,"_compare_",aa,"_VS_",bb,".pdf"), width = 10,height = 8.5)
    print(gmat)
    dev.off()
  }
}





message("donut")