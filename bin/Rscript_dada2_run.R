library(seqinr)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(tibble)
library(foreach)
suppressMessages(library(doParallel))

ncores<-as.numeric(commandArgs()[7])
subp<-commandArgs()[8]
fwdname<-commandArgs()[9]
rvsname<-commandArgs()[10]
minov<-as.numeric(commandArgs()[11])
maxmis<-minov*(1-as.numeric(commandArgs()[12]))
prev<-commandArgs()[13]

# all libraries
librun<-read.table("../config/librun.list",sep="\t",col.names=c("sample","library","run"),stringsAsFactors=F)

# input libraries
lib<-ddply(librun,.(run,library,sample),function(x){
  fwdin<-list.files(path=x$sample,pattern=paste0(x$library,".fwd.filtered.fastq"),full.names=T)
  rvsin<-list.files(path=x$sample,pattern=paste0(x$library,".rvs.filtered.fastq"),full.names=T)
  rbind(data.frame(lib="fwd",filename=fwdin,stringsAsFactors=F),
        data.frame(lib="rvs",filename=rvsin,stringsAsFactors=F))
}) %>%
  mutate(dir=sub("\\..*$","",basename(filename)))

# priors ASV
if (prev != "no") {
  prior_asv<-data.frame(filename=list.files(pattern=paste0(prev,".*\\.[fr][wv][ds]\\.fasta$")),stringsAsFactors=F) %>%
    separate(filename,c("subp","dir","lib","ext"),sep="\\.",remove=F)
  prior_seq<-dlply(prior_asv,.(lib),function(x) {
    unique(foreach(y=iapply(x,1),.combine=c) %do% {
      unlist(read.fasta(y$filename,seqonly=T))
    })
  })
}

# error model
err<-dlply(lib,.(run,lib), function(x) dada2::learnErrors(x$filename,randomize=T,multithread=ncores))

# dada2
ncpus<-ifelse(nrow(lib)<ncores,nrow(lib),ncores)
cl<-makeCluster(ncpus)
registerDoParallel(cl)
dada_all<-dlply(lib, .(run,lib),function(x) {
  tmp_err<-err[[unique(paste(x$run,x$lib,sep="."))]]
  if (prev != "no") {
    tmp_prior<-prior_seq[[unique(x$lib)]]
    dlply(x,.(sample,library,dir),function(y) {
      tmp_derep<-dada2::derepFastq(y$filename,n=1e5)
      dada2::dada(tmp_derep,tmp_err,priors=tmp_prior)
    },.parallel=T)
  } else {
    dlply(x,.(sample,library,dir),function(y) {
      tmp_derep<-dada2::derepFastq(y$filename,n=1e5)
      dada2::dada(tmp_derep,tmp_err)
    },.parallel=T)
  }
})
stopCluster(cl)

# merge pairs
pairs<-mutate(lib,comb=ifelse((lib=="fwd" & dir==fwdname) | (lib=="rvs" & dir==rvsname),"FR","RF")) %>%
  pivot_wider(names_from=lib,values_from=c(dir,filename))
dada_pairs<-dlply(pairs,.(run,library,sample,comb),function(x) {
  list(sample=x$sample,
       fwd=dada_all[[paste0(x$run,".fwd")]][[paste(x$sample,x$library,x$dir_fwd,sep=".")]],
       filename_fwd=x$filename_fwd,
       rvs=dada_all[[paste0(x$run,".rvs")]][[paste(x$sample,x$library,x$dir_rvs,sep=".")]],
       filename_rvs=x$filename_rvs)
})
rm(dada_all)
cl<-makeCluster(ncores)
registerDoParallel(cl)
mergers_all<-llply(dada_pairs,function(x) {
  tmp_derep_F<-dada2::derepFastq(file.path(x$filename_fwd),n=1e5)
  tmp_derep_R<-dada2::derepFastq(file.path(x$filename_rvs),n=1e5)
  dada2::mergePairs(x$fwd, tmp_derep_F, x$rvs, tmp_derep_R, minOverlap=minov, maxMismatch=maxmis, trimOverhang=T)
},.parallel=T,.paropts=list(.export=c('minov','maxmis')))
stopCluster(cl)

# sequence table for each paired library in each sample, sum per sample, reverse-coomplement if RF direction
seqtab_pairs<-llply(setNames(as.list(unique(pairs$comb)),unique(pairs$comb)),function(x) {
  dada2::mergeSequenceTables(tables=llply(setNames(as.list(unique(pairs$sample)),unique(pairs$sample)),function(y) {
    tmp<-colSums(dada2::makeSequenceTable(mergers_all[grep(paste0(y,"\\.",x,"$"),names(mergers_all),value=T)]))
    if(x=="RF") {
      return(matrix(tmp,nrow=1,dimnames=list(paste(y,x,sep="_"),
                                             laply(names(tmp),function(x) c2s(rev(comp(s2c(x),forceToLower=F)))))))
    } else {
      return(matrix(tmp,nrow=1,dimnames=list(paste(y,x,sep="_"),names(tmp))))
    }
  }))
})

# merge to a single table
if (length(seqtab_pairs)>1) {
  seqtab<-dada2::mergeSequenceTables(tables=seqtab_pairs)
} else {
  seqtab<-seqtab_pairs[[1]]
}

# Remove chimeras
seqtab_clean<-dada2::removeBimeraDenovo(seqtab, method="pooled", multithread=ncores)

# percent of ASV and reads removed
nbbim<-ncol(seqtab)-ncol(seqtab_clean)
print(paste0("Found ",nbbim," chimera (",round(nbbim/ncol(seqtab)*100,digits=2)," %) and ",
             ncol(seqtab_clean)," (",round(ncol(seqtab_clean)/ncol(seqtab)*100,digits=2)," %) non-chimera ASVs."))
countbim<-sum(seqtab)-sum(seqtab_clean)
print(paste0("Found ",countbim," chimera (",round(countbim/sum(seqtab)*100,digits=2)," %) and ",
             sum(seqtab_clean)," (",round(sum(seqtab_clean)/sum(seqtab)*100,digits=2)," %) non-chimera amplicons."))

# merge sequence counts from both sequencing direction
mat_df<-setNames(data.frame(seqtab_clean),paste0("ASV_",sprintf(paste0("%0",nchar(ncol(seqtab_clean)),"d"),1:ncol(seqtab_clean)))) %>%
  rownames_to_column("sample") %>%
  mutate(sample=sub("_[FR][FR]$","",sample)) %>%
  group_by(sample) %>%
  summarize_all(sum)
final_asv<-data.frame(asv=colnames(mat_df)[-1],seq=colnames(seqtab_clean),stringsAsFactors=F)

# transpose
mat<-column_to_rownames(mat_df,"sample") %>%
  t() %>%
  data.frame(total=colSums(mat_df[,-1]),.) %>%
  rownames_to_column("Representative_Sequence")

# rename ASV after previous ASV, if any
if (prev != "no") {
  prev_asv<-read.fasta(file.path(paste0(prev,".outputs"),paste0(prev,".all_repseq.fasta")),as.string=T,forceDNAtolower=F,set.attributes=F) %>%
    ldply(function(x) c(seq=x),.id="prev") %>%
    mutate(prev=as.character(prev)) %>%
    inner_join(final_asv,.,by="seq") %>%
    select(-seq) %>%
    arrange(prev)
  prev_length<-length(grep("^>",readLines(file.path(paste0(prev,".outputs"),paste0(prev,".all_repseq.fasta")))))
  mat<-mutate(mat,asv=Representative_Sequence) %>%
    right_join(prev_asv,.,by="asv") %>%
    arrange(prev,-total) %>%
    separate(asv,c("header","index"),sep="_") %>%
    mutate(index=as.numeric(index),
           Representative_Sequence=ifelse(is.na(prev),
                                          paste0("ASV_",sprintf(paste0("%0",nchar(n()+prev_length),"d"),index+prev_length)),
                                          prev)) %>%
    select(-header,-index,-prev)
  final_asv<-rbind(right_join(final_asv,prev_asv,by="asv") %>% transmute(asv=prev,seq=seq),
                   filter(final_asv,! asv %in% prev_asv$asv) %>%
                     separate(asv,c("header","index"),sep="_",remove=F) %>%
                     mutate(index=as.numeric(index),
                            asv=paste0("ASV_",sprintf(paste0("%0",nchar(n()+prev_length),"d"),index+prev_length))) %>%
                     select(-header,-index))
}

# track all sequences
map_track<-dlply(pairs,.(run,library,sample,comb),function(x) {
  tmp_derep_F<-dada2::derepFastq(file.path(x$filename_fwd),n=1e5)
  tmp_derep_R<-dada2::derepFastq(file.path(x$filename_rvs),n=1e5)
  tmp_dada_F<-dada_pairs[[paste(x$run,x$library,x$sample,x$comb,sep=".")]][["fwd"]]
  tmp_dada_R<-dada_pairs[[paste(x$run,x$library,x$sample,x$comb,sep=".")]][["rvs"]]
  tmp_mergers<-mergers_all[[paste(x$run,x$library,x$sample,x$comb,sep=".")]]
  tmp_map<-full_join(data.frame(filtered=1:length(tmp_derep_F$map),derep_forward=tmp_derep_F$map,derep_reverse=tmp_derep_R$map),
                     data.frame(derep_forward=1:length(tmp_dada_F$map),dada_forward=tmp_dada_F$map),by="derep_forward") %>%
    full_join(data.frame(derep_reverse=1:length(tmp_dada_R$map),dada_reverse=tmp_dada_R$map),by="derep_reverse") %>%
    full_join(data.frame(dada_forward=tmp_mergers$forward,dada_reverse=tmp_mergers$reverse,merged=1:length(tmp_mergers$forward),
                         seq=if(x$comb=="RF"){laply(tmp_mergers$sequence,function(y) c2s(rev(comp(s2c(y),forceToLower=F))))} else {tmp_mergers$sequence},
                         stringsAsFactors=F),by=c("dada_forward","dada_reverse")) %>%
    full_join(final_asv,by="seq") %>%
    filter(!is.na(filtered)) %>%
    select(-seq)
})
map_track_sample<-ldply(map_track) %>%
  group_by(run,library,sample,comb) %>%
  summarise_all(~length(na.exclude(.))) %>%
  group_by(sample) %>%
  select(-run,-library,-comb,-starts_with("derep")) %>%
  summarize_all(~sum(.)) %>%
  rename(bimera_removed=asv)

# export ASVs as fasta, the ASV table with named ASVs, the track index for each filtered sequence to each ASV and count statistic
write(apply(final_asv,1,function(x) paste0(">",paste(x,collapse="\n"))),file=paste(subp,"dada2.fasta",sep="."),ncolumns=1)
write.table(mat,paste0(subp,".dada2.count_table"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(map_track_sample,paste0(subp,".dada2.read_counts.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
it<-as.list(iapply(pairs,1))
for(i in 1:length(it)) {
  x<-it[[i]]
  write.table(select(map_track[[paste(x$run,x$library,x$sample,x$comb,sep=".")]],filtered,asv) %>%
                filter(!is.na(asv)),
              file.path(sub("\\.fastq$",".asv.index",x$filename_fwd)),col.names=F,row.names=F,quote=F)
}

