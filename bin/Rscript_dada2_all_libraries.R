library(seqinr)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(tibble)
library(foreach)
suppressMessages(library(doParallel))
suppressMessages(library(dada2))

ncores<-as.numeric(commandArgs()[7])
subp<-commandArgs()[8]
fwdname<-commandArgs()[9]
rvsname<-commandArgs()[10]
minov<-as.numeric(commandArgs()[11])
maxmis<-minov*(1-as.numeric(commandArgs()[12]))
prev<-commandArgs()[13]

# libraries and samples
lib3<-read.table("../config/lib3.list",sep="\t",stringsAsFactors=F)
samples<-unique(lib3[,1])
libraries<-foreach(i=samples,.combine=rbind) %do% {
  foreach(j=c("fwd","rvs"),.combine=rbind) %do% {
    foreach(k=c(fwdname,rvsname),.combine=rbind) %do% {
      tmp<-list.files(path=i,pattern=paste(k,".*",j,"filtered.derep",sep="."))
      if (length(tmp)>0) {
        data.frame(sample=i,lib=j,dir=k,filename=tmp,stringsAsFactors=F) %>%
          mutate(library=sub(paste0("^",k,"\\."),"",sub(paste0("\\.",j,"\\.filtered\\.derep$"),"",filename)))
      }
    }
  }
}
pairs<-mutate(libraries,comb=ifelse((lib=="fwd" & dir==fwdname) | (lib=="rvs" & dir==rvsname),"FR","RF")) %>%
  pivot_wider(names_from=lib,values_from=c(dir,filename))

# priors ASV
prior_asv<-data.frame(filename=list.files(pattern="*.[fr][wv][ds].fasta$"),stringsAsFactors=F) %>%
  separate(filename,c("dir","lib","ext"),sep="\\.",remove=F)
prior_seq<-dlply(prior_asv,.(lib),function(x) {
  foreach(y=iapply(x,1),.combine=c) %do% {
    unlist(read.fasta(y$filename,seqonly=T))
  }
})

# dada pseudo-pooling stategy with prior ASV from all libraries and track sequence map
cl<-makeCluster(ncores)
registerDoParallel(cl)
dada_pairs<-suppressWarnings(dlply(pairs,.(sample,library,comb),function(x) {
  tmp_err_fwd<-readRDS(file.path(x$sample,paste0(x$library,".fwd.err.rds")))
  tmp_derep_fwd<-readRDS(file.path(x$sample,x$filename_fwd))
  tmp_err_rvs<-readRDS(file.path(x$sample,paste0(x$library,".rvs.err.rds")))
  tmp_derep_rvs<-readRDS(file.path(x$sample,x$filename_rvs))
  list(sample=x$sample,
       fwd=dada(tmp_derep_fwd,tmp_err_fwd,priors=prior_seq$fwd),
       filename_fwd=x$filename_fwd,
       rvs=dada(tmp_derep_rvs,tmp_err_rvs,priors=prior_seq$fwd),
       filename_rvs=x$filename_rvs)
},.parallel=T,.paropts=list(.export='prior_seq',.packages='dada2')))
stopCluster(cl)

# merge pairs
cl<-makeCluster(ncores)
registerDoParallel(cl)
mergers_all<-suppressWarnings(llply(dada_pairs,function(x) {
  tmp_derep_fwd<-readRDS(file.path(x$sample,x$filename_fwd))
  tmp_derep_rvs<-readRDS(file.path(x$sample,x$filename_rvs))
  mergePairs(x$fwd, tmp_derep_fwd, x$rvs, tmp_derep_rvs, minOverlap=minov, maxMismatch=maxmis, trimOverhang=T)
},.parallel=T,.paropts=list(.export=c('minov','maxmis'),.packages='dada2')))
stopCluster(cl)

# sequence table for each paired library in each sample, sum per sample, reverse-coomplement if RF direction
seqtab_pairs<-llply(setNames(as.list(unique(pairs$comb)),unique(pairs$comb)),function(x) {
  mergeSequenceTables(tables=llply(setNames(as.list(samples),samples),function(y) {
    tmp<-colSums(makeSequenceTable(mergers_all[grep(paste(y,".*",x,sep="\\."),names(mergers_all),value=T)]))
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
  seqtab<-mergeSequenceTables(tables=seqtab_pairs)
} else {
  seqtab<-seqtab_pairs[[1]]
}

# Remove chimeras
seqtab_clean<-removeBimeraDenovo(seqtab, method="pooled", multithread=ncores)

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
  data.frame(total=rowSums(.),.,check.names=F) %>%
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
map_track<-dlply(pairs,.(sample,library,comb),function(x) {
  tmp_derep_F<-readRDS(file.path(x$sample,x$filename_fwd))
  tmp_derep_R<-readRDS(file.path(x$sample,x$filename_rvs))
  tmp_dada_F<-dada_pairs[[paste(x$sample,x$library,x$comb,sep=".")]][["fwd"]]
  tmp_dada_R<-dada_pairs[[paste(x$sample,x$library,x$comb,sep=".")]][["rvs"]]
  tmp_mergers<-mergers_all[[paste(x$sample,x$library,x$comb,sep=".")]]
  tmp_map<-full_join(data.frame(libraries=1:length(tmp_derep_F$map),derep_forward=tmp_derep_F$map,derep_reverse=tmp_derep_R$map),
                            data.frame(derep_forward=1:length(tmp_dada_F$map),dada_forward=tmp_dada_F$map),by="derep_forward") %>%
    full_join(data.frame(derep_reverse=1:length(tmp_dada_R$map),dada_reverse=tmp_dada_R$map),by="derep_reverse") %>%
    full_join(data.frame(dada_forward=tmp_mergers$forward,dada_reverse=tmp_mergers$reverse,merged=1:length(tmp_mergers$forward),
                         seq=if(x$comb=="RF"){laply(tmp_mergers$sequence,function(y) c2s(rev(comp(s2c(y),forceToLower=F))))} else {tmp_mergers$sequence},
                         stringsAsFactors=F),by=c("dada_forward","dada_reverse")) %>%
    full_join(final_asv,by="seq") %>%
    filter(!is.na(libraries)) %>%
    select(-seq)
})
map_track_sample<-ldply(map_track) %>%
  group_by(sample,library,comb) %>%
  summarise_all(~length(na.exclude(.))) %>%
  group_by(sample) %>%
  select(-library,-comb,-starts_with("derep"),-seq) %>%
  summarize_all(~sum(.)) %>%
  rename(bimera_removed=asv)

# export ASVs as fasta, the ASV table with named ASVs, the track index for each libraries sequence to each ASV and count statistic
write(apply(final_asv,1,function(x) paste0(">",paste(x,collapse="\n"))),file=paste(subp,"dada2.fasta",sep="."),ncolumns=1)
write.table(mat,paste0(subp,".dada2.count_table"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(map_track_sample,paste0(subp,".dada2.read_counts.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
it<-as.list(iapply(pairs,1))
for(i in 1:length(it)) {
  x<-it[[i]]
  write.table(select(map_track[[paste(x$sample,x$library,x$comb,sep=".")]],libraries,asv) %>%
    filter(!is.na(asv)),
    file.path(x$sample,sub("\\.derep$",".asv.index",x$filename_fwd)),col.names=F,row.names=F,quote=F)
}
