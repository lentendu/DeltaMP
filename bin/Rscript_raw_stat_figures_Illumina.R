library(plyr)
suppressMessages(library(dplyr))
library(tidyr)
library(tibble)
library(foreach)
library(ggplot2)
subp<-commandArgs()[7]
bin<-commandArgs()[8]
amp<-commandArgs()[9]
lib<-commandArgs()[10]
fwd<-commandArgs()[11]
rvs<-commandArgs()[12]
trunc<-commandArgs()[13]
source(file.path(bin,"palette.R"))

if(trunc=="no") {
	types<-c("fwd","rvs","pairend","pairend_trim")
	fulltypes<-c("forward reads","reverse reads","pair-end reads","minimal quality paired reads")
} else {
	types<-c("fwd","rvs","pairend")
	fulltypes<-c("forward reads","reverse reads","pair-end reads")
}
list<-apply(expand.grid(types,c("meanqual","meanposqual","length")),1,paste,collapse=".")

if(lib=="all") {
  assign(paste0("files_",fwd),sub(paste0(fwd,"\\."),"",sub("\\.fwd\\.meanqual","",list.files(pattern=paste0(fwd,".*.fwd.meanqual")))))
  assign(paste0("files_",rvs),sub(paste0(rvs,"\\."),"",sub("\\.rvs\\.meanqual","",list.files(pattern=paste0(rvs,".*.rvs.meanqual")))))
  samples<-unique(sort(sub("\\..*$","",c(get(paste0("files_",fwd)),get(paste0("files_",rvs))))))
  lname<-subp
  project<-paste0("Project: ",subp)
} else {
  assign(paste0("files_",fwd),sub(paste0(fwd,"\\."),"",sub("\\.fwd\\.meanqual","",list.files(pattern=paste0(fwd,".*",lib,".fwd.meanqual")))))
  assign(paste0("files_",rvs),sub(paste0(rvs,"\\."),"",sub("\\.rvs\\.meanqual","",list.files(pattern=paste0(rvs,".*",lib,".rvs.meanqual")))))
  samples<-unique(c(sub(paste0("\\.",lib),"",get(paste0("files_",fwd))),sub(paste0("\\.",lib),"",get(paste0("files_",rvs)))))
  lname<-paste(subp,lib,sep=".")
  project<-paste0("Project: ",subp,"\nLibrary: ",lib)
}

mycolors<-palette()[1:length(samples)]
my_theme<-theme_bw() +
  theme(strip.background=element_rect(fill="grey90",linetype=0),
        strip.text=element_text(size=12),
        panel.border=element_blank(),
        axis.line=element_line(color="black"),
        plot.title=element_text(size=16,hjust=0.5),
        plot.subtitle=element_text(size=12,hjust=0.5),
        plot.margin=margin(t=12, r=12, b=12, l=12, unit="pt"))
theme_set(my_theme)

meanqual<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    foreach(z=c(fwd,rvs),w=c(rvs,fwd),.combine=rbind) %do% {
      if(x %in% c("fwd","rvs")) {
        tmp<-paste(z,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"meanqual",sep=".")
      } else {
        tmp<-paste(z,w,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"meanqual",sep=".")
      }
      if(file.exists(tmp)) {
        data.frame(type=x,sample=y,dir=z,quality=0:41,reads=scan(tmp,quiet=T))
      }
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes))

meanposqual<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    foreach(z=c(fwd,rvs),w=c(rvs,fwd),.combine=rbind) %do% {
      if(x %in% c("fwd","rvs")) {
        tmp<-paste(z,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"meanposqual",sep=".")
      } else {
        tmp<-paste(z,w,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"meanposqual",sep=".")
      }
      if(file.exists(tmp)) {
        data.frame(type=x,sample=y,dir=z,quality=scan(tmp,quiet=T)) %>%
          rownames_to_column("position") 
      }
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes),position=as.numeric(position))

seqlength<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    foreach(z=c(fwd,rvs),w=c(rvs,fwd),.combine=rbind) %do% {
      if(x %in% c("fwd","rvs")) {
        tmp<-paste(z,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"length",sep=".")
      } else {
        tmp<-paste(z,w,grep(paste0("^",y),get(paste0("files_",z)),value=T),x,"length",sep=".")
      }
      if(file.exists(tmp)) {
        data.frame(type=x,sample=y,dir=z,read.table(tmp,col.names=c("length","reads"))) 
      }
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes)) %>%
  filter(length>0)

# average quality
for (x in c(fwd,rvs)) {
	if (length(get(paste0("files_",x))) > 0) {
	  plot_qual<-ggplot(filter(meanqual,dir==x),aes(quality,reads+1)) +
	    geom_bar(data=group_by(filter(meanqual,dir==x),type,quality) %>% summarize(reads=sum(reads)),stat="identity",fill="limegreen") +
	    geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
	    geom_smooth() +
	    facet_wrap(~type,ncol=1) +
	    scale_color_manual(values=mycolors) +
	    scale_y_log10() +
	    labs(title=paste0("Raw reads quality for reads with ",x," primer at 5'-end"),subtitle=paste0(project,"\n",amp,paste(rep(" ",58),collapse=" "),Sys.Date()),
	         x="quality score",y="read counts")
	  ggsave(paste(lname,"quality",x,"pdf",sep="."),plot_qual,width=210,height=297,units="mm")
	  
	  # average quality by position
	  plot_pos_qual<-ggplot(filter(meanposqual,dir==x),aes(position,quality)) +
	    geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
	    geom_smooth() +
	    facet_wrap(~type,ncol=1,scales="free_x") +
	    scale_color_manual(values=mycolors) +
	    labs(title=paste0("Length vs. quality distribution for reads with ",x," primer at 5'-end"),x="nucleotide position",y="quality score")
	  ggsave(paste(lname,"position_quality",x,"pdf",sep="."),plot_pos_qual,width=210,height=297,units="mm")
	  
	  # sequence length
	  plot_length<-ggplot(filter(seqlength,dir==x),aes(length,reads)) +
	    geom_bar(data=group_by(filter(seqlength,dir==x),length,type) %>%
	               summarize(reads=sum(reads)),stat="identity",fill="limegreen") +
	    geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
	    geom_smooth() +
	    facet_wrap(~type,ncol=1,scales="free_x") +
	    scale_color_manual(values=mycolors) +
	    scale_y_log10() +
	    labs(title=paste0("length distribution for reads with ",x," primer at 5'-end"),x="sequence length (nt)",y="read counts")
	  ggsave(paste(lname,"length",x,"pdf",sep="."),plot_length,width=210,height=297,units="mm")
  }
}
# legend
plot_line<-ggplot(seqlength,aes(length,reads)) +
  geom_line(aes(color=sample),stat="identity") +
  scale_color_manual(values=mycolors,guide=guide_legend(ncol=3)) +
  theme(legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.spacing.x=unit(5, 'mm'),
        legend.key.size=unit(5,'mm'))
plot_legend<-cowplot::get_legend(plot_line)
pdf(paste(lname,"legend","pdf",sep="."),paper="a4",width=0,height=0)
plot(plot_legend)
dev.off()
