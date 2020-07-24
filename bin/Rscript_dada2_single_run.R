library(seqinr)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(tibble)
library(foreach)
suppressMessages(library(doParallel))

ncores<-as.numeric(commandArgs()[7]) # ncores<-3
runname<-commandArgs()[8] # runname<-"test_lib1.R"
fwdname<-commandArgs()[9] # fwdname<-"TAReuk454FWD1"
rvsname<-commandArgs()[10] # rvsname<-"TAReukREV3"
comb<-commandArgs()[11] # comb<-"FR"
minov<-as.numeric(commandArgs()[12]) # minov<-10
maxmis<-minov*(1-as.numeric(commandArgs()[13])) # maxmis<-minov*(1-0.6)
prev<-commandArgs()[14] # prev<-"no"

# all libraries
librun<-read.table("../config/librun.list",sep="\t",col.names=c("sample","library","run"),stringsAsFactors=F) %>%
  filter(run==runname) %>%
  select(-run)

# input libraries
libs<-ddply(librun,.(sample),function(x){
  fwdin<-list.files(path=x$sample,pattern=paste(fwdname,x$library,"fwd.filtered.fastq",sep="."),full.names=T)
  rvsin<-list.files(path=x$sample,pattern=paste(rvsname,x$library,"rvs.filtered.fastq",sep="."),full.names=T)
  rbind(data.frame(lib="fwd",filename=fwdin,stringsAsFactors=F),
        data.frame(lib="rvs",filename=rvsin,stringsAsFactors=F))
}) %>%
  mutate(dir=sub("\\..*$","",basename(filename)))

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

# error model for all samples of the same run and and same lib (R1/R2)
err<-dlply(libs,.(lib),function(x) {
  set.seed(1)
  dada2::learnErrors(x$filename,randomize=T,multithread=ncores)
})

# dada2: pool strategy for all samples of same run, same lib and same orientation (i.e. starting by forward or reverse primer)
ncpus<-ifelse(nrow(libs)<ncores,nrow(libs),ncores)
dada_all<-dlply(libs,.(lib),function(x) {
  tmp_err<-err[[unique(x$lib)]]
  cl<-makeCluster(ncpus)
  registerDoParallel(cl)
  tmp_derep<-suppressWarnings(plyr::dlply(x,.(sample),function(y) {
    dada2::derepFastq(y$filename,n=1e5)
  },.parallel=T))
  stopCluster(cl)
 if (prev != "no") {
   tmp_prior<-prior_seq[[unique(paste(y$dir,y$lib,sep="."))]]
   dada2::dada(tmp_derep[grep(unique(y$dir),names(tmp_derep),value=T)],tmp_err,pool=T,multithread=ncores,priors=tmp_prior)
 } else {
   dada2::dada(tmp_derep,tmp_err,pool=T,multithread=ncores)
 }
})

# merge pairs
pairs<-pivot_wider(libs,names_from=lib,values_from=c(dir,filename))
dada_pairs<-dlply(pairs,.(sample),function(x) {
  list(sample=x$sample,
       fwd=dada_all[["fwd"]][[x$sample]],
       filename_fwd=x$filename_fwd,
       rvs=dada_all[["rvs"]][[x$sample]],
       filename_rvs=x$filename_rvs)
})
rm(dada_all)
gc()
cl<-makeCluster(ncores)
registerDoParallel(cl)
mergers_all<-suppressWarnings(llply(dada_pairs,function(x) {
  tmp_derep_F<-dada2::derepFastq(file.path(x$filename_fwd),n=1e5)
  tmp_derep_R<-dada2::derepFastq(file.path(x$filename_rvs),n=1e5)
  dada2::mergePairs(x$fwd, tmp_derep_F, x$rvs, tmp_derep_R, minOverlap=minov, maxMismatch=maxmis, trimOverhang=T)
},.parallel=T,.paropts=list(.export=c('minov','maxmis'))))
stopCluster(cl)

# sequence table for each paired library in each sample, sum per sample, reverse-complement if RF direction
seqtab_pairs<-dada2::makeSequenceTable(mergers_all)
if(comb=="RF") {
  colnames(seqtab_pairs)<-laply(colnames(seqtab_pairs),function(x) c2s(rev(comp(s2c(x),forceToLower=F))))
}
rownames(seqtab_pairs)<-paste(rownames(seqtab_pairs),runname,comb,sep="_")

# track all sequences
map_track<-dlply(pairs,.(sample),function(x) {
  tmp_derep_F<-dada2::derepFastq(file.path(x$filename_fwd),n=1e5)
  tmp_derep_R<-dada2::derepFastq(file.path(x$filename_rvs),n=1e5)
  tmp_dada_F<-dada_pairs[[x$sample]][["fwd"]]
  tmp_dada_R<-dada_pairs[[x$sample]][["rvs"]]
  tmp_mergers<-mergers_all[[x$sample]]
  tmp_map<-full_join(data.frame(filtered=1:length(tmp_derep_F$map),derep_forward=tmp_derep_F$map,derep_reverse=tmp_derep_R$map),
                     data.frame(derep_forward=1:length(tmp_dada_F$map),dada_forward=tmp_dada_F$map),by="derep_forward") %>%
    full_join(data.frame(derep_reverse=1:length(tmp_dada_R$map),dada_reverse=tmp_dada_R$map),by="derep_reverse") %>%
    full_join(data.frame(dada_forward=tmp_mergers$forward,dada_reverse=tmp_mergers$reverse,merged=1:length(tmp_mergers$forward),
                         seq=if(comb=="RF"){laply(tmp_mergers$sequence,function(y) c2s(rev(comp(s2c(y),forceToLower=F))))} else {tmp_mergers$sequence},
                         stringsAsFactors=F),by=c("dada_forward","dada_reverse"))
})
map_track_sample<-ldply(map_track) %>%
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
  write.table(select(map_track[[x$sample]],filtered,seq) %>%
                filter(!is.na(seq)),
              file.path(sub("\\.fastq$",".asv.index",x$filename_fwd)),col.names=F,row.names=F,quote=F)
}

