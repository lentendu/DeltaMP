library(seqinr)
library(digest)
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dada2))

ncores<-as.numeric(commandArgs()[7])
subp<-commandArgs()[8]
fwdname<-commandArgs()[9]
rvsname<-commandArgs()[10]
minov<-as.numeric(commandArgs()[11])
maxmis<-minov*(1-as.numeric(commandArgs()[12]))
minlen<-as.numeric(commandArgs()[13])
bim<-commandArgs()[14]
pseudo<-commandArgs()[15]

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

cl<-makeCluster(ncores)
registerDoParallel(cl)
if ( pseudo == "yes" ) {
	# priors ASV
	prior_asv<-data.frame(filename=list.files(pattern="*.[fr][wv][ds].fasta$"),stringsAsFactors=F) %>%
	  separate(filename,c("dir","lib","ext"),sep="\\.",remove=F)
	prior_seq<-dlply(prior_asv,.(lib),function(x) {
	  foreach(y=iapply(x,1),.combine=c) %do% {
	    unlist(read.fasta(y$filename,seqonly=T))
	  }
	})
	
	# dada pseudo-pooling stategy with prior ASV from all libraries and track sequence map
	dada_pairs<-suppressWarnings(dlply(pairs,.(sample,library,comb),function(x) {
	  tmp_err_fwd<-readRDS(file.path(x$sample,paste0(x$library,".fwd.err.rds")))
	  tmp_derep_fwd<-readRDS(file.path(x$sample,x$filename_fwd))
	  tmp_err_rvs<-readRDS(file.path(x$sample,paste0(x$library,".rvs.err.rds")))
	  tmp_derep_rvs<-readRDS(file.path(x$sample,x$filename_rvs))
	  list(sample=x$sample,
	       fwd=dada2::dada(tmp_derep_fwd,tmp_err_fwd,priors=prior_seq$fwd),
	       filename_fwd=x$filename_fwd,
	       rvs=dada2::dada(tmp_derep_rvs,tmp_err_rvs,priors=prior_seq$fwd),
	       filename_rvs=x$filename_rvs)
	},.parallel=T,.paropts=list(.export='prior_seq')))
} else {
  dada_pairs<-suppressWarnings(dlply(pairs,.(sample,library,comb),function(x) {
    list(sample=x$sample,
         fwd=readRDS(file.path(x$sample,sub("filtered.derep","dada",x$filename_fwd))),
         filename_fwd=x$filename_fwd,
         rvs=readRDS(file.path(x$sample,sub("filtered.derep","dada",x$filename_rvs))),
         filename_rvs=x$filename_rvs)
  },.parallel=T))
}

# merge pairs
mergers_all<-suppressWarnings(llply(dada_pairs,function(x) {
  tmp_derep_fwd<-readRDS(file.path(x$sample,x$filename_fwd))
  tmp_derep_rvs<-readRDS(file.path(x$sample,x$filename_rvs))
  dada2::mergePairs(x$fwd, tmp_derep_fwd, x$rvs, tmp_derep_rvs, minOverlap=minov, maxMismatch=maxmis, trimOverhang=T)
},.parallel=T,.paropts=list(.export=c('minov','maxmis'))))
stopCluster(cl)

# list samples with sequences in each orientation (avoid processing orientations without sequences)
sample_pairs<-dlply(attr(dada_pairs,"split_labels"),.(comb))

# sequence table for each paired library in each sample, sum per sample, reverse-coomplement if RF direction
seqtab_pairs<-llply(setNames(as.list(unique(pairs$comb)),unique(pairs$comb)),function(x) {
  mergeSequenceTables(tables=llply(getElement(sample_pairs,x)$sample %>% as.list() %>% setNames(.,.),function(y) {
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

# Remove too short sequences
seqtab_ok<-seqtab[,nchar(colnames(seqtab))>=minlen]

# Remove bimeras
if ( bim == "yes" ) {
	seqtab_clean<-removeBimeraDenovo(seqtab_ok, method="consensus", multithread=ncores)
} else {
	seqtab_clean<-seqtab_ok
}

# percent of ASV and reads removed
nbbim<-ncol(seqtab)-ncol(seqtab_clean)
print(paste0("Found ",nbbim," bimera (",round(nbbim/ncol(seqtab)*100,digits=2)," %) and ",
             ncol(seqtab_clean)," (",round(ncol(seqtab_clean)/ncol(seqtab)*100,digits=2)," %) non-bimera ASVs."))
countbim<-sum(seqtab)-sum(seqtab_clean)
print(paste0("Found ",countbim," bimera (",round(countbim/sum(seqtab)*100,digits=2)," %) and ",
             sum(seqtab_clean)," (",round(sum(seqtab_clean)/sum(seqtab)*100,digits=2)," %) non-bimera amplicons."))

# merge sequence counts from both sequencing direction
cl<-makeCluster(ncores)
registerDoParallel(cl)
mat_df<-setNames(data.frame(seqtab_clean),laply(colnames(seqtab_clean),sha1)) %>%
  rownames_to_column("sample") %>%
  mutate(sample=sub("(_.*)*_[FR][FR]$","",sample)) %>%
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
    filter(!is.na(libraries))
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
cl<-makeCluster(ncores)
registerDoParallel(cl)
idout<-foreach(x=iapply(pairs,1), .packages=c('dplyr')) %dopar% {
  write.table(filter(map_track[[paste(x$sample,x$library,x$comb,sep=".")]],!is.na(seq)) %>% rowwise() %>% mutate(sha1=digest::sha1(seq)) %>% select(libraries,sha1),
    file.path(x$sample,sub("\\.derep$",".asv.index",x$filename_fwd)),col.names=F,row.names=F,quote=F)
}
stopCluster(cl)
