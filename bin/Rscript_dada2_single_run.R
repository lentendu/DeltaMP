library(seqinr)
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dada2))

ncores<-as.numeric(commandArgs()[7])
runname<-commandArgs()[8]
fwdname<-commandArgs()[9]
rvsname<-commandArgs()[10]
comb<-commandArgs()[11]
minov<-as.numeric(commandArgs()[12])
maxmis<-minov*(1-as.numeric(commandArgs()[13]))
prev<-commandArgs()[14]

# all libraries
librun<-read.table("../config/librun.list",sep="\t",col.names=c("sample","library","run"),stringsAsFactors=F) %>%
  filter(run==runname) %>%
  select(-run)

# input libraries
libs<-ddply(librun,.(sample),function(x){
  fwdin<-list.files(path=x$sample,pattern=paste(fwdname,x$library,"fwd.filtered.derep",sep="."),full.names=T)
  rvsin<-list.files(path=x$sample,pattern=paste(rvsname,x$library,"rvs.filtered.derep",sep="."),full.names=T)
  rbind(data.frame(lib="fwd",filename=fwdin,stringsAsFactors=F),
        data.frame(lib="rvs",filename=rvsin,stringsAsFactors=F))
}) %>%
  mutate(dir=sub("\\..*$","",basename(filename)))
pairs<-pivot_wider(libs,names_from=lib,values_from=c(dir,filename))


# priors ASV
if (prev != "no") {
  prior_asv<-data.frame(filename=list.files(pattern=paste0(prev,".*\\.[fr][wv][ds]\\.fasta$")),stringsAsFactors=F) %>%
    separate(filename,c("subp","dir","lib","ext"),sep="\\.",remove=F)
  prior_seq<-dlply(prior_asv,.(dir,lib),function(x) {
    unique(foreach(y=iapply(x,1),.combine=c) %do% {
      unlist(read.fasta(y$filename,seqonly=T))
    })
  })
}

# dereplicated
derep<-dlply(libs,.(lib),function(x) dlply(x,.(sample), function(y) readRDS(y$filename)))

# error model for all samples of the same run and same lib (R1/R2)
err<-dlply(libs,.(lib),function(x) {
  learnErrors(derep[[unique(x$lib)]],multithread=ncores)
})

# dada2: pool strategy for all samples of same run, same lib and same orientation (i.e. starting by forward or reverse primer)
dada_all<-dlply(libs,.(lib),function(x) {
  tmp_err<-err[[unique(x$lib)]]
  tmp_derep<-derep[[unique(x$lib)]]
 if (prev != "no") {
   tmp_prior<-prior_seq[[unique(paste(y$dir,y$lib,sep="."))]]
   dada(tmp_derep,tmp_err,pool=T,multithread=ncores,priors=tmp_prior)
 } else {
   dada(tmp_derep,tmp_err,pool=T,multithread=ncores)
 }
})

# merge pairs
cl<-makeCluster(ncores)
registerDoParallel(cl)
mergers_all<-foreach(i=dada_all$fwd,j=derep$fwd,k=dada_all$rvs,l=derep$rvs,.packages='dada2',
                     .final=function(x)setNames(x,names(derep$fwd))) %dopar% {
                       mergePairs(i,j,k,l, minOverlap=minov, maxMismatch=maxmis, trimOverhang=T)
                     }
stopCluster(cl)

# sequence table for each paired library in each sample, sum per sample, reverse-complement if RF direction
seqtab_pairs<-makeSequenceTable(mergers_all)
if(comb=="RF") {
  colnames(seqtab_pairs)<-laply(colnames(seqtab_pairs),function(x) c2s(rev(comp(s2c(x),forceToLower=F))))
}
rownames(seqtab_pairs)<-paste(rownames(seqtab_pairs),runname,comb,sep="_")

# track all sequences
cl<-makeCluster(ncores)
registerDoParallel(cl)
map_track<-foreach(i=dada_all$fwd,j=derep$fwd,k=dada_all$rvs,l=derep$rvs,m=mergers_all,.packages=c('plyr','dplyr','seqinr'),
                     .final=function(x)setNames(x,names(derep$fwd))) %dopar% {
                       full_join(data.frame(filtered=1:length(j$map),derep_forward=j$map,derep_reverse=l$map),
                                 data.frame(derep_forward=1:length(i$map),dada_forward=i$map),by="derep_forward") %>%
                         full_join(data.frame(derep_reverse=1:length(k$map),dada_reverse=k$map),by="derep_reverse") %>%
                         full_join(data.frame(dada_forward=m$forward,dada_reverse=m$reverse,merged=1:length(m$forward),
                                              seq=if(comb=="RF"){
                                                laply(m$sequence,function(y) c2s(rev(comp(s2c(y),forceToLower=F))))
                                                } else {m$sequence},
                                              stringsAsFactors=F),by=c("dada_forward","dada_reverse"))
                     }
stopCluster(cl)
map_track_sample<-ldply(map_track,.id="sample") %>%
  select(-seq) %>%
  group_by(sample) %>%
  summarise_all(~length(na.exclude(.))) %>%
  select(-starts_with("derep")) %>%
  mutate(run=runname,comb=comb)

# export
saveRDS(seqtab_pairs,paste(runname,comb,"dada2_seqtab.rds",sep="##"))
saveRDS(map_track_sample,paste(runname,comb,"dada2_track.rds",sep="##"))
it<-as.list(iapply(pairs,1))
for(i in 1:length(it)) {
  x<-it[[i]]
  write.table(select(map_track[[x$sample]],filtered,seq) %>% filter(!is.na(seq)),
              file.path(sub("\\.derep$",".asv.index",x$filename_fwd)),col.names=F,row.names=F,quote=F)
}

