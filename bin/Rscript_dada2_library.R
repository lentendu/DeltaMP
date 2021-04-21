library(plyr)
library(foreach)
suppressMessages(library(doParallel))
suppressMessages(library(dada2))

ncores<-as.numeric(commandArgs()[7])
lib<-commandArgs()[8]
perrun<-commandArgs()[9]

# input libraries
fwdin<-list.files(pattern=paste0(lib,".fwd.filtered.fastq"))
rvsin<-list.files(pattern=paste0(lib,".rvs.filtered.fastq"))

# avoid empty files
for (i in c("fwdin","rvsin")) {
  assign(i,get(i)[file.info(get(i))$size>0])
}

# dereplicate and export
frnames<-rbind(cbind.data.frame(lib="fwd",dir=sub(paste0("\\.",lib,".*$"),"",fwdin),stringsAsFactors=F),
               cbind.data.frame(lib="rvs",dir=sub(paste0("\\.",lib,".*$"),"",rvsin),stringsAsFactors=F))
if(nrow(frnames)>ncores) {ncpus<-ncores} else {ncpus<-nrow(frnames)}
cl<-makeCluster(ncpus)
registerDoParallel(cl)
derep<-suppressWarnings(dlply(frnames, .(lib,dir),function(i) {
  tmpf<-grep(i$dir,get(paste0(i$lib,"in")),value=T)
  tmp<-derepFastq(tmpf,n=1e5)
  saveRDS(tmp,sub("fastq$","derep",tmpf))
  return(tmp)
},.parallel=T,.paropts=list(.export=c('fwdin','rvsin'),.packages='dada2')))
stopCluster(cl)


# start pseudo-pooling routine if not per run
if ( perrun == "no" ) {
  # error model
  err<-dlply(frnames,.(lib),function(i) {
    set.seed(1)
    learnErrors(derep[apply(i,1,paste,collapse=".")],randomize=T,multithread=ncores)
  })
  
  # dada2
  dada<-dlply(frnames, .(lib,dir),
              function(i) dada(derep[[paste(i$lib,i$dir,sep=".")]],err[[i$lib]],multithread=ncores))
  
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
}

