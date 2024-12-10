suppressMessages(library(plyr))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

nhits<-commandArgs()[7]
cons<-as.numeric(commandArgs()[8])
ncores<-as.numeric(commandArgs()[9])
srange<-as.numeric(commandArgs()[10])

hits<-read.table(nhits,sep="\t",col.names=c("seqid","sim","tax")) %>%
  separate(tax,c("ref","taxo"),sep=";tax=",fill="right") %>%
  mutate(taxo=sub(";$","",taxo)) %>% 
  arrange(seqid,desc(sim)) %>%
  filter(sim >= max(sim)-((100-max(sim))*srange),.by=seqid)

cl<-makeCluster(ncores)
registerDoParallel(cl)
suppressWarnings(
  out<-ddply(hits,.(seqid),function(y) {
    print(unique(y$seqid))
    if ( all(is.na(y$taxo)) ) {
      data.frame(seqid=y$seqid,sim=0,taxo="no_match",ref=NA)
    } else {
      tmp<-separate(y,taxo,paste0("taxo_",c(letters,LETTERS)),sep=",",fill="right") %>%
        select_if(~!(all(is.na(.)))) %>%
        select(starts_with("taxo_"))
      ntaxo=NULL ; nboot=NULL
      for(i in 1:ncol(tmp)) {
        perc<-round(table(tmp[i])/nrow(tmp)*100)
        if(any(perc>cons)) {
          ntaxo[i]<-names(perc[which(perc>cons)])
          nboot[i]<-perc[ntaxo[i]]  
        } else {
          break 
        }
      }
      if (is.null(ntaxo)) {
        data.frame(seqid=unique(y$seqid),
                   sim=NA,
                   taxo=NA,
                   ref=NA)
      } else {
        tmp2<-filter(y,grepl(ntaxo[length(ntaxo)],taxo))
        data.frame(seqid=unique(y$seqid),
                   sim=ifelse(length(unique(tmp2$sim))>1,paste(min(tmp2$sim),max(tmp2$sim),sep="-"),unique(tmp2$sim)),
                   taxo=sub("$",";",paste(paste0(ntaxo,"(",nboot,")"),collapse=";")),
                   ref=paste(tmp2$ref,collapse=","))
      }
    }
  },.parallel=T,.paropts=list(.packages=c('tidyr','dplyr','foreach'),.export=c('cons')))
)
stopCluster(cl)

write.table(out,sub("hits$","taxonomy",nhits),sep="\t",row.names=F,col.names=F,quote=F)
