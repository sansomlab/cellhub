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

swapR <- function(matrixRow,x,y){
  #x is diagonal index
  #y is max of the row
  indexY <- which(matrixRow == y)
  valX <- matrixRow[x]
  matrixRow[x] <- y
  matrixRow[indexY] <- valX
  return(matrixRow)
}




option_list <- list(
  make_option(c("--basedir"), default="/well/singlecell/P190381/_zzz_/10X-count",
              help="the folder where project data of the demultiplexing. the name has to finish with _10x-count-_ to match with Santiago's new structure"),
  make_option(c("--demultiplexing"), default="demuxlet,vireo,demuxlet2",
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

if (!file.exists(opt$outdir)){ dir.create(opt$outdir)}
print("Running with options:")
print(opt)

#parse adn format outdir, print opts statement
run<- opt$outdir
opt[which(opt=="NULL")] <- NULL


opt$demultiplexing<-strsplit(opt$demultiplexing,",")[[1]]
gsub(" ","",opt$demultiplexing)->opt$demultiplexing
if(length(opt$demultiplexing)<2){ #uniform names results
  if (opt$demultiplexing=="demuxlet" | opt$demultiplexing=="demuxlet2"){ 
    dend<- ".best"
  }else if(opt$demultiplexing=="vireo") {
    dend<- "donor_ids.tsv"
  }else{
    stop("can't recognise this multiplexing output")
  }
  pstring <- paste0( opt$basedir,"/", opt$demultiplexing, "/",opt$samplename,"*/*",dend)
  ll <- system(paste0("ls ",pstring), intern=T)
  dlist <-list()
  if(opt$demultiplexing=="demuxlet") {
    cid<-c("BARCODE","BEST")
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      assign(paste0("demuxlet",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("demuxlet",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
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
  
  if(opt$demultiplexing=="demuxlet2") {
    cid<-c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      assign(paste0("demuxletV2",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("demuxlet",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
      general<-dlist[[idx]][,2]
      specific<-dlist[[idx]][,3]
      which(general=="SNG") ->nex
      general[nex]<-specific[nex]
      which(general == "AMB") ->NEX
      general[NEX] <- "UNASSIGNED"
      which(general == "DBL") ->DEX
      general[DEX] <- "DOUBLET"
      dlist[[idx]][,2] <- general
      nex <- NEX <-DEX <- NULL
      dlist[[idx]]<-dlist[[idx]][,c(1,2)]
    }
    names(dlist) <- basename(dirname(ll)) 
  }
  
  
  
  if(opt$demultiplexing=="vireo"){
    cid<-c("cell","donor_id")
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      assign(paste0("vireo",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("vireo",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
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
  if("demuxlet"  %in%  opt$demultiplexing){
  dm<-opt$demultiplexing[which(opt$demultiplexing=="demuxlet")]
  dstring<- paste0( opt$basedir,"/demuxlet/",opt$samplename,"*/*.best")
  }
  if("demuxlet2"  %in%  opt$demultiplexing){
  dm2<-opt$demultiplexing[which(opt$demultiplexing=="demuxlet")]
  d2string<- paste0( opt$basedir,"/demuxlet2/",opt$samplename,"*/*.best")
  }
  if("vireo"  %in%  opt$demultiplexing){
  vm<-opt$demultiplexing[which(opt$demultiplexing=="vireo")]
  vstring<- paste0( opt$basedir, "/vireo/",opt$samplename,"*/donor_ids.tsv")
  }
  
  if("demuxlet"  %in%  opt$demultiplexing){
    message("demuxlet")
    demtemp <- "demuxlet"
    cid<-c("BARCODE","BEST")
    ll <- system(paste0("ls ",dstring), intern=T)
    dlist <-list()
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      assign(paste0("demuxlet",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      colnames(dlist[[idx]])[2] <- paste0("demuxlet",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
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
  message("Data demultiplexed with Demuxlet v1:")
  print(dim(data.demuxlet))
  }
  
  
  if("demuxlet2"  %in%  opt$demultiplexing){
    message("demuxletV2")
    cid<-c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")
    ll <- system(paste0("ls ",d2string), intern=T)
    dlist <-dem2temp<-list()
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      assign(paste0("demuxletV2",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      dem2temp[[idx]] <- dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("demuxletV2",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
      colnames(dem2temp[[idx]])[2] <- paste0("demuxletV2",gsub(opt$samplename,"",basename(dirname(ll[idx]))) )
      colnames(dem2temp[[idx]])[3] <- paste0("demuxletV2",gsub(opt$samplename,"",basename(dirname(ll[idx]))),"SNG.BEST.GUESS") 
      
      general<-dlist[[idx]][,2]
      specific<-dlist[[idx]][,3]
      which(general=="SNG") ->nex
      general[nex]<-specific[nex]
      which(general == "AMB") ->NEX
      general[NEX] <- "UNASSIGNED"
      which(general == "DBL") ->DEX
      general[DEX] <- "DOUBLET"
      dlist[[idx]][,2] <- general
      dlist[[idx]]<-dlist[[idx]][,c(1,2)]
      nex <- NEX <-DEX <- NULL
    }
    names(dem2temp) <-basename(dirname(ll)) 
    names(dlist) <- basename(dirname(ll)) 
    data.demuxletv2<- Reduce(function(x,y) merge(x,y,by="BARCODE"), dlist)
    message("Data demultiplexed with Demuxlet v2:")
    print(dim(data.demuxletv2))
    #---------------------------------
    data.temp.demuxletv2<-Reduce(function(x,y) merge(x,y,by="BARCODE"), dem2temp)
    output_path <-  file.path(paste0(run, opt$samplename, "allmethods.demuxlet2.best.doublet.tsv"))
    write.table(data.temp.demuxletv2,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)
    gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)
    for( NN in 1:length(dem2temp)){
      ppst<-names(dem2temp[[NN]])
      output_path <-  file.path(paste0(run, ppst,"_demuxlet2.best.doublet.tsv"))
      write.table(dem2temp[[NN]],file =output_path, sep="\t", quote = F, row.names = F, col.names = T)
    gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)
    
    }
    rm( data.temp.demuxletv2,dem2temp,NN,pst)
  }
  


  
  if("vireo"  %in%  opt$demultiplexing){
    message("vireo")
    demtemp <- "vireo"
    cid<-c("cell","donor_id")
    ll <- system(paste0("ls ",vstring), intern=T)
    dlist <-list()
    for( idx in 1:length(ll)){
      if (idx %% 2 == 0) {message("file ", idx, " of ", length(ll), " ...")}
      assign(paste0("vireo",gsub(opt$samplename,"",basename(dirname(ll[idx])))),read.table(ll[idx], header=T))
      dlist[[idx]]<-read.table(ll[idx], header=T)[,cid]
      colnames(dlist[[idx]])[2] <- paste0("vireo",gsub(opt$samplename,"",basename(dirname(ll[idx])))) 
      colnames(dlist[[idx]])[1] <-"BARCODE"
      dlist[[idx]][,2]<-toupper(dlist[[idx]][,2])
    }
    names(dlist) <- basename(dirname(ll)) 
    data.vireo<- Reduce(function(x,y) merge(x,y,by="BARCODE"), dlist)
    message("Data demultiplexed with vireo:")
    print(dim(data.vireo))
    
  }
  
  xtt<-list() 
  for (N in 1:length(ls()[grepl("data.", ls())])) { xtt[[N]] <-get (ls()[grepl("data.", ls())][N]) }
  
  data.use <-Reduce(function(x,y) merge(x,y,by="BARCODE"), xtt) 
  
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

infomat %>% 
  pivot_wider(id_cols=value,
              values_from = tot, 
              names_from = variable) %>% 
  replace(is.na(.),0) ->save.sng.cnt



infomat %>% mutate(General=ifelse(!(value %in% c("UNASSIGNED","DOUBLET")), "SINGLET",value )) %>% 
  group_by(variable, General) %>% summarise(tot=sum(tot)) %>%
  group_by(variable) %>% mutate(persample=sum(tot), percent=100*tot/persample) -> infomat.general

infomat.general %>% 
  pivot_wider(id_cols=General,
              values_from = percent, 
              names_from = variable) %>% 
  replace(is.na(.),0) ->save.general

output_path <-  file.path(paste0(run, opt$samplename, "_SingleCellMetadata_demultiplexing_results.tsv"))
write.table(data.use,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)

output_path <-  file.path(paste0(run, opt$samplename, "_General_demultiplexing_results.tsv"))
write.table(save.general,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)

output_path <-  file.path(paste0(run, opt$samplename, "_Singlets_demultiplexing_results.tsv"))
write.table(save.sng,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)


output_path <-  file.path(paste0(run, opt$samplename, "_Singlets_demultiplexing_results_cellcount.tsv"))
write.table(save.sng.cnt,file =output_path, sep="\t", quote = F, row.names = T, col.names = T)
gzip(output_path,destname=sprintf("%s.gz", output_path), overwrite=TRUE, remove=TRUE)


c("#7d00b3ff","#93ff00ff","#4e79e9ff") ->scol


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
c("DOUBLET","UNASSIGNED",gnr[!gnr %in% c("DOUBLET","UNASSIGNED")])->gnr
sum(!gnr %in% c("DOUBLET","UNASSIGNED")) ->vc
colormap(colormap=colormaps$portland,vc) ->ccol
c("#7d00b3ff","#93ff00ff", ccol) ->vcol


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

infomat %>% mutate(value=factor(value, levels = gnr)) %>% 
  ggplot(aes(variable, tot, fill=value)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values = vcol, 
                    breaks=gnr, 
                    limits=gnr) +
  theme(axis.text.x = element_text(angle=72, hjust=1 ,vjust=1),
        legend.position="top", legend.key.size = unit(0.3, "cm")) +
  guides(colour = guide_legend(nrow = 3)) + ylab("cell yield") ->g1

ggsave(g1, filename = file.path(paste0(run, opt$samplename, "_Singlets_demultiplexing_results_cellCount.pdf")), width=10, height = 7)



colnames(data.use)[!colnames(data.use)=="BARCODE"]->cols

ddf<-matrix(0,nrow=length(cols), ncol = length(cols))
rownames(ddf)<- colnames(ddf)<- cols
col_fun = colorRamp2(c(0, 0.50, 1.00), c("dodgerblue4", 'peachpuff', 'deeppink4'))

for(a in cols){
  for(b in cols){
    ddf[a,b] <- mclust::adjustedRandIndex(data.use[,a],data.use[,b])  
  }
}
  
Heatmap(ddf, cluster_rows = F,cluster_columns = F, col = col_fun,
        row_names_side = "left", row_names_gp = gpar(fontsize=9),
        column_names_gp = gpar(fontsize=9),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.4f", ddf[i, j]), x, y, gp = gpar(fontsize = 8))
        }) ->gg

pdf(paste0(run, opt$samplename,"_adjustedRand_scores.pdf"), 10,10)
gg
dev.off()


ddf %>% as.data.frame() %>% rownames_to_column("Var1") %>% 
  pivot_longer(cols, names_to = "Var2") %>% 
  filter(Var2 >Var1) %>% 
  arrange(-value) %>% rename(AdjustedRandIndex=value)->arranged.rand

output_path <-  file.path(paste0(run, opt$samplename, "_RandIndexScores.txt"))
write.table(arranged.rand,file =output_path, sep="\t", quote = F, row.names = F, col.names = T)


ddf<-matrix(0,nrow=length(cols), ncol = length(cols))
rownames(ddf)<- colnames(ddf)<- cols

unique(unlist(data.use[,cols])) ->allpeop
allpeop[!allpeop  %in% c("DOUBLET","UNASSIGNED")] ->ALL



for(a in cols){
  for(b in cols){
    
    data.use %>% select(a,b) %>% 
      filter_all(all_vars(. !="DOUBLET")) %>% 
      filter_all(all_vars(. !="UNASSIGNED")) ->snt
    mat.all <-matrix(0, ncol = length(ALL), nrow=length(ALL))
    rownames(mat.all) <-colnames(mat.all) <- ALL
    
    mat<-table(snt[,a],snt[,b])
    mat.all[rownames(mat.all) %in% rownames(mat), 
            colnames(mat.all) %in% colnames(mat)] <- mat[rownames(mat) %in% rownames(mat.all), 
                                                         colnames(mat) %in% colnames(mat.all)]
    
    mat.max <- matrix(0, nrow=nrow(mat.all), ncol=ncol(mat.all))
    for(i in 1:nrow(mat.all)){
      rowI <- mat.all[i,]
      y <- max(rowI)
      mat.max[i,] <- swapR(rowI, i, y)
    }
    
    ddf[a,b] <- sum(diag(  mat.max))/sum(  mat.max)
  }
}

col_fun = colorRamp2(c(0.50,0.80, 1.00), c("dodgerblue4", 'peachpuff', 'deeppink4'))
Heatmap(ddf, cluster_rows = F,cluster_columns = F, col = col_fun,
        row_names_side = "left", row_names_gp = gpar(fontsize=9),
        column_names_gp = gpar(fontsize=9),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.4f", ddf[i, j]), x, y, gp = gpar(fontsize = 8))
        }) ->gg

pdf(paste0(run, opt$samplename,"_SingletConcordance.pdf"), 10,10)
gg
dev.off()

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
  gsub(" ", "",opt$subset) ->opt$subset
  data.use %>% select(ends_with(c(opt$subset, as.character(best$variable) ))) ->data.subset
}else {
  data.use %>% select(ends_with(as.character(best$variable) )) ->data.subset
}

col_fun = colorRamp2(c(0, 50, 100), c("dodgerblue4", 'peachpuff', 'deeppink4'))
dir.create(paste0(run,"comparison_selection"))
cols <- colnames(data.subset)
for(a in cols){
  for(b in cols){
    if (a !=b){
    n<-colSums(table(data.subset[,a],data.subset[,b])) #colsums over B
    norm.mat<-t(apply(table(data.subset[,a],data.subset[,b]),1, function(x) 100*x/n))
    Heatmap(norm.mat, 
            cluster_rows = F,
            cluster_columns = F,
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
            name = paste0("Percent\nover\n",b),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", norm.mat[i, j]), x, y, gp = gpar(fontsize = 7))
            }) ->gmat
    aa <- gsub(paste0(".",opt$samplename),"",a)
    bb <- gsub(paste0(".",opt$samplename),"",b)
    pdf(paste0(run,"comparison_selection/",opt$samplename,"_compare_",aa,"_VS_",bb,".pdf"), width = 10,height = 8.5)
    print(gmat)
    dev.off()
    }
  }
}


colnames(data.subset) ->vartest
message("plotting selected samples")
for(demtype in vartest){
  message(demtype)
  if (grepl("vireo",demtype)){
    get(demtype) %>% 
      filter(cell  %in% data.use$BARCODE) %>% 
      ggplot(aes(donor_id, n_vars, color=prob_doublet, fill=donor_id)) +
      geom_violin() + geom_jitter(size=1) + theme_linedraw() + 
      theme(axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12,angle=72, hjust=1, vjust=1   ),
            legend.position = "bottom", legend.key.size = unit(0.5, "cm")) +
      scale_fill_manual(values = vcol, 
                        breaks= toupper(gnr), 
                        limits= toupper(gnr)) +
      scale_color_viridis_c(option="A") +
      ylab("Number of Variants")->g
    ggsave(g, filename = file.path(paste0(run, opt$samplename, demtype,"demreport.pdf")), width=11, height = 8)    
      
  }else if (grepl("^demuxlet$",demtype)){
    general<-sapply(strsplit(get(demtype)[,2],"-"),"[[",1)
    specific<-sapply(strsplit(get(demtype)[,2],"-"),"[[",2)
    which(general=="SNG") ->nex
    general[nex]<-specific[nex]
    which(general == "AMB") ->NEX
    general[NEX] <- "UNASSIGNED"
    which(general == "DBL") ->DEX
    general[DEX] <- "DOUBLET"
    get(demtype)->temp
    
    temp[,"general"]<-general
    assign(demtype, temp)
    nex <- NEX <-DEX <-temp<- NULL
    
    get(demtype) %>% 
      filter(BARCODE  %in% data.use$BARCODE) %>% 
      ggplot(aes(general, N.SNP, color=PRB.DBL, fill=general)) +
      geom_violin() + geom_jitter(size=1) + theme_linedraw() + 
      theme(axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12,angle=72, hjust=1, vjust=1),
            legend.position = "bottom", legend.key.size = unit(0.5, "cm"))+
      scale_fill_manual(values = vcol, 
                        breaks= toupper(gnr), 
                        limits= toupper(gnr)) +
      scale_color_viridis_c(option="A") +
      ylab("Number of Variants")->g
    ggsave(g, filename = file.path(paste0(run, opt$samplename, demtype,"demreport.pdf")), width=11, height = 8)    
      
  }else if (grepl("^demuxletV2",demtype)){
    
    general<-get(demtype)[,"DROPLET.TYPE"]
    specific<-get(demtype)[,"SNG.BEST.GUESS"]
    which(general=="SNG") ->nex
    general[nex]<-specific[nex]
    which(general == "AMB") ->NEX
    general[NEX] <- "UNASSIGNED"
    which(general == "DBL") ->DEX
    general[DEX] <- "DOUBLET"
    get(demtype)->temp
    
    temp[,"general"]<-general
    assign(demtype, temp)
    nex <- NEX <-DEX <-temp<- NULL
    
    
    
    get(demtype) %>%  
      filter(BARCODE  %in% data.use$BARCODE) %>% 
      ggplot(aes(general, NUM.SNPS, fill=general)) +
      geom_violin() + geom_jitter(size=1) + theme_linedraw() + 
      theme(axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12,angle=72, hjust=1, vjust=1),
            legend.position = "bottom", legend.key.size = unit(0.5, "cm"))+
      scale_fill_manual(values = vcol, 
                        breaks= toupper(gnr), 
                        limits= toupper(gnr)) +
      
      ylab("Number of Variants")->g
    ggsave(g, filename = file.path(paste0(run, opt$samplename, demtype,"demreport.pdf")), width=11, height = 8)    
    
    
  }
}




message("donut")