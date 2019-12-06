library(plyr)
suppressMessages(library(dplyr))
library(tidyr)
library(tibble)
library(foreach)
library(ggplot2)
suppressMessages(library(gridExtra))
subp<-commandArgs()[7]
bin<-commandArgs()[8]
amp<-commandArgs()[9]
lib<-commandArgs()[10]
source(file.path(bin,"palette.R"))

list<-apply(expand.grid(c("fwd","rvs","pairend","pairend_trim"),c("meanqual","meanposqual","length")),1,paste,collapse=".")

if(lib=="all") {
  files<-sub(paste0("\\.fwd\\.meanqual"),"",list.files(pattern="fwd.meanqual"))
  samples<-files
  lname<-subp
  project<-paste0("Project: ",subp)
} else {
  files<-sub(paste0("\\.fwd\\.meanqual"),"",list.files(pattern=paste(lib,"fwd.meanqual",sep=".")))
  samples<-sub(paste0("\\.",lib),"",files)
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

types<-c("fwd","rvs","pairend","pairend_trim")
fulltypes<-c("forward reads","reverse reads","pair-end reads","minimal quality paired reads")

meanqual<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    tmp<-paste(grep(paste0("^",y,"$"),files,value=T),x,"meanqual",sep=".")
    if(file.exists(tmp)) {
      data.frame(type=x,sample=y,quality=0:41,reads=scan(tmp,quiet=T))
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes))

meanposqual<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    tmp<-paste(grep(paste0("^",y,"$"),files,value=T),x,"meanposqual",sep=".")
    if(file.exists(tmp)) {
      data.frame(type=x,sample=y,quality=scan(tmp,quiet=T)) %>%
        rownames_to_column("position") 
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes),position=as.numeric(position))

seqlength<-foreach(x=types,.combine=rbind) %do% {
  foreach(y=samples,.combine=rbind) %do% {
    tmp<-paste(grep(paste0("^",y,"$"),files,value=T),x,"length",sep=".")
    if(file.exists(tmp)) {
      data.frame(type=x,sample=y,read.table(tmp,col.names=c("length","reads")))
    }
  }
} %>% mutate(type=mapvalues(type,types,fulltypes)) %>%
  filter(length>0)

# average quality
plot_qual<-ggplot(meanqual,aes(quality,reads+1)) +
  geom_bar(data=group_by(meanqual,type,quality) %>% summarize(reads=sum(reads)),stat="identity",fill="limegreen") +
  geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
  geom_smooth() +
  facet_wrap(~type,ncol=1) +
  scale_color_manual(values=mycolors) +
  scale_y_log10() +
  labs(title="Raw reads quality",subtitle=paste0(project,"\n",amp,paste(rep(" ",58),collapse=" "),Sys.Date()),
       x="quality score",y="read counts")
ggsave(paste(lname,"quality","pdf",sep="."),plot_qual,width=210,height=297,units="mm")

# average quality by position
plot_pos_qual<-ggplot(meanposqual,aes(position,quality)) +
  geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
  geom_smooth() +
  facet_wrap(~type,ncol=1,scales="free_x") +
  scale_color_manual(values=mycolors) +
  labs(title="Length vs. quality distribution",x="nucleotide position",y="quality score")
ggsave(paste(lname,"position_quality","pdf",sep="."),plot_pos_qual,width=210,height=297,units="mm")


# sequence length
plot_length<-ggplot(seqlength,aes(length,reads)) +
  geom_bar(data=group_by(seqlength,length,type) %>%
             summarize(reads=sum(reads)),stat="identity",fill="limegreen") +
  geom_line(aes(color=sample),size=0.3,stat="identity",show.legend=F) +
  geom_smooth() +
  facet_wrap(~type,ncol=1,scales="free_x") +
  scale_color_manual(values=mycolors) +
  scale_y_log10() +
  labs(title="length distribution",x="sequence length (nt)",y="read counts")
ggsave(paste(lname,"length","pdf",sep="."),plot_length,width=210,height=297,units="mm")

# legend
plot_line<-ggplot(seqlength,aes(length,reads)) +
  geom_line(aes(color=sample),stat="identity") +
  scale_color_manual(values=mycolors,guide=guide_legend(ncol=2)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.spacing.x=unit(5, 'mm'),
        legend.key.size=unit(5,'mm'))
plot_legend<-cowplot::get_legend(plot_line)
pdf(paste(lname,"legend","pdf",sep="."),paper="a4",width=0,height=0)
plot(plot_legend)
dev.off()
