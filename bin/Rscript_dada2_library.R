library(plyr)
library(foreach)
suppressMessages(library(doParallel))

ncores<-as.numeric(commandArgs()[7])
lib<-commandArgs()[8]

# input libraries
fwdin<-list.files(pattern=paste0(lib,".fwd.filtered.fastq"))
rvsin<-list.files(pattern=paste0(lib,".rvs.filtered.fastq"))

# error model
err<-llply(setNames(as.list(c("fwd","rvs")),c("fwd","rvs")),function(i) {
  set.seed(1)
  dada2::learnErrors(get(paste0(i,"in")),randomize=T,multithread=ncores)
})

# dereplicate
frnames<-rbind(cbind.data.frame(lib="fwd",dir=sub(paste0("\\.",lib,".*$"),"",fwdin),stringsAsFactors=F),
               cbind.data.frame(lib="rvs",dir=sub(paste0("\\.",lib,".*$"),"",rvsin),stringsAsFactors=F))
cl<-makeCluster(nrow(frnames))
registerDoParallel(cl)
derep<-dlply(frnames, .(lib,dir),
             function(i) dada2::derepFastq(grep(i$dir,get(paste0(i$lib,"in")),value=T),n=1e5),
             .parallel=T,.paropts=list(.export=c('fwdin','rvsin')))
stopCluster(cl)

# dada2
dada<-dlply(frnames, .(lib,dir),
            function(i) dada2::dada(derep[[paste(i$lib,i$dir,sep=".")]],err[[i$lib]],multithread=ncores))

# export error models and ASVs
for (i in c("fwd","rvs")) {
  saveRDS(err[[i]],paste(lib,i,"err.rds",sep="."))
  for (j in unique(frnames$dir)) {
    write(paste0(">",paste(j,lib,i,1:length(dada[[paste(i,j,sep=".")]]$sequence),sep="."),
                 "\n",dada[[paste(i,j,sep=".")]]$sequence),
          paste(j,lib,i,"asv.fasta",sep="."),
          ncolumns=1) 
  }
}

