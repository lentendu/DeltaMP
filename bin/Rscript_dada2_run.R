library(seqinr)
library(digest)
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dada2))

ncores<-as.numeric(commandArgs()[7])
subp<-commandArgs()[8]
minlen<-as.numeric(commandArgs()[9])

# merge all seqtab and remove too short sequences
list_seqtab<-data.frame(seqtab=list.files(pattern=".dada2_seqtab.rds$"),stringsAsFactors=F) %>%
  separate(seqtab,c("lib","comb","suf"),sep="##",remove=F) %>%
  select(-suf)
seqtab<-dlply(list_seqtab,.(lib,comb),function(x) readRDS(x$seqtab) %>% .[,nchar(colnames(.))>=minlen]) %>%
  mergeSequenceTables(tables=.)

# remove bimeras
seqtab_clean<-removeBimeraDenovo(seqtab, method="pooled", multithread=ncores)

# percent of ASV and reads removed
nbbim<-ncol(seqtab)-ncol(seqtab_clean)
print(paste0("Found ",nbbim," chimera (",round(nbbim/ncol(seqtab)*100,digits=2)," %) and ",
             ncol(seqtab_clean)," (",round(ncol(seqtab_clean)/ncol(seqtab)*100,digits=2)," %) non-chimera ASVs."))
countbim<-sum(seqtab)-sum(seqtab_clean)
print(paste0("Found ",countbim," chimera (",round(countbim/sum(seqtab)*100,digits=2)," %) and ",
             sum(seqtab_clean)," (",round(sum(seqtab_clean)/sum(seqtab)*100,digits=2)," %) non-chimera amplicons."))

# merge sequence counts from both sequencing direction (and from multiple runs)
cl<-makeCluster(ncores)
registerDoParallel(cl)
mat_df<-setNames(data.frame(seqtab_clean),laply(colnames(seqtab_clean),sha1)) %>%
  rownames_to_column("sample") %>%
  separate(sample,c("sample","run","dir"),sep="@") %>%
  select(-run,-dir) %>%
  ddply(.(sample),function(x){
    if(is.null(dim(x))) {x} else {colSums(dplyr::select(x,-sample))}
  },.parallel=T)
stopCluster(cl)
final_asv<-data.frame(asv=colnames(mat_df)[-1],seq=colnames(seqtab_clean),stringsAsFactors=F)

# transpose
mat<-column_to_rownames(mat_df,"sample") %>%
  t() %>%
  data.frame(total=rowSums(.),.,check.names=F) %>%
  rownames_to_column("Representative_Sequence")

# stat on sequences
list_track<-data.frame(track=list.files(pattern=".dada2_track.rds$"),stringsAsFactors=F) %>%
  separate(track,c("xlib","xcomb","suf"),sep="##",remove=F) %>%
  select(-suf)
map_track_sample<-ddply(list_track,.(xlib,xcomb),function(x) readRDS(x$track)) %>%
  select(-xlib,-xcomb,-run,-comb) %>%
  group_by(sample) %>%
  summarize_all(~sum(.)) %>%
  full_join(data.frame(bimera_removed=colSums(mat[,-c(1,2)])) %>%
              rownames_to_column("sample"))

# export ASVs as fasta, the ASV table with named ASVs, the track index for each filtered sequence to each ASV and count statistic
write(apply(final_asv,1,function(x) paste0(">",paste(x,collapse="\n"))),file=paste(subp,"dada2.fasta",sep="."),ncolumns=1)
write.table(mat,paste0(subp,".dada2.count_table"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(map_track_sample,paste0(subp,".dada2.read_counts.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
