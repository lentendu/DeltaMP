library(plyr)
subp<-commandArgs()[7]
bin<-commandArgs()[8]
amp<-commandArgs()[9]
source(file.path(bin,"palette.R"))
 
# Read stat tables
list1<-sub("\\.stat","",list.files(pattern="*\\.stat$"))
for (i in list1) {
  assign(i, read.table(paste0(i,".stat"),h=T,check.names=F))
}

# Add position/length column if necessary
maxc<-max(aaply(list1,1,function(x) ncol(get(x))))
for (i in list1) {
  if(ncol(get(i))<maxc) {
    tmp<-cbind.data.frame(seq(1,nrow(get(i))),get(i))
    colnames(tmp)[1]<-i
    assign(i,tmp)
  }
}

# sum
for (i in list1) {
  tmp<-cbind.data.frame(get(i)[,1],sum=round(apply(get(i)[,-1],1,sum)))
  colnames(tmp)<-c(colnames(get(i))[1],"sum")
  assign(paste0("sum.",i),tmp)
}

# Figures

## Raw reads ##
pdf(paste(subp,"raw_reads_with_primer_quality.pdf",sep="."),paper="a4",width=0,height=0,title=paste(subp,"raw reads quality"))
layout(matrix(c(1:5),ncol=1),heights=c(1,4,4,4,4))
# Title
par(mar=c(0,2,0,0),cex=1)
plot.new()
text(0.5,0.85,"Raw reads quality (with primers)",font=2,xpd=T,cex=1.4)
text(0,0.3,paste0("Project: ",subp),adj=c(0,NA),font=2,xpd=T,cex=1)
text(1,0.3,paste0(amp,"   ",Sys.Date()),adj=c(1,NA),xpd=T,cex=1)
## raw forward and reverse ##
labfr<-c(fwd="forward",rvs="reverse")
# Mean sequence quality
par(mar=c(3,3,1,0.5),mgp=c(1.8,0.6,0),tck=-0.05,cex.axis=0.8,las=1,xpd=T)
for (i in c("fwd","rvs")) {
  tmp_qual<-get(paste0(i,".meanqual"))
  tmp_sum<-get(paste0("sum.",i,".meanqual"))
  lmax<-max(tmp_sum[,2])
  posx<-unique(floor(fwd.meanqual[,1]/5))*5
  posy<-10^(1:nchar(lmax))
  barplot(log(tmp_sum[,2]+1),axes=F,border=NA,space=0,ylim=c(0,log(posy[length(posy)]+1)),col="limegreen",xlab=paste(labfr[i],"average phred score"))
  axis(1,at=seq(posx[1]-tmp_qual[1,1]+1,which(tmp_qual[,1]==posx[length(posx)]),5),labels=posx)
  axis(2,at=log(c(0,posy)+1),labels=c(0,sub("\\+0","\\+",format(posy,scientific=T))))
  mtext("Count",side=2,xpd=T,adj=0.5,padj=-3.7,las=0)
  for (i in 2:ncol(tmp_qual)) {
    points(seq(0.5,nrow(tmp_qual)-0.5,1),log(tmp_qual[,i]+1),type="l",col=palette()[i]) 
  }
}
# Mean base quality per position
for (i in c("fwd","rvs")) {
  tmp_qual<-get(paste0(i,".meanposqual"))
  plot(NULL,bty="n",xlim=c(0,nrow(tmp_qual)),ylim=c(floor(min(tmp_qual[,-1])/10)*10,40),xlab=paste(labfr[i],"nucleotide position"),ylab="Phred score")
  for (i in 2:ncol(tmp_qual)) {
    points(tmp_qual[,1],tmp_qual[,i],type="l",col=palette()[i]) 
  }
}
dev.off()


## Pair-end ##
pdf(paste(subp,"pair-end_reads_quality.pdf",sep="."),paper="a4",width=0,height=0,title=paste(subp,"pair-end reads quality"))
layout(matrix(c(1:5),ncol=1),heights=c(1,4,4,4,4))
# Title
par(mar=c(0,2,0,0),cex=1)
plot.new()
text(0.5,0.85,"Pair-end reads quality",font=2,xpd=T,cex=1.4)
text(0,0.3,paste0("Project: ",subp),adj=c(0,NA),font=2,xpd=T,cex=1)
text(1,0.3,paste0(amp,"   ",Sys.Date()),adj=c(1,NA),xpd=T,cex=1)
# Mean sequence quality
par(mar=c(3,3,1,0.5),mgp=c(1.8,0.5,0),tck=-0.03,cex.axis=0.8,las=1,xpd=T)
lmax<-max(sum.pairend.meanqual[,2])
posx<-unique(floor(fwd.meanqual[,1]/5))*5
posy<-10^(1:nchar(lmax))
barplot(log(sum.pairend.meanqual[,2]+1),axes=F,border=NA,space=0,ylim=c(0,log(posy[length(posy)]+1)),col="limegreen",xlab="Average phred score")
axis(1,at=seq(posx[1]-pairend.meanqual[1,1]+1,which(pairend.meanqual[,1]==posx[length(posx)]),5),labels=posx)
axis(2,at=log(c(0,posy)+1),labels=c(0,sub("\\+0","\\+",format(posy,scientific=T))))
mtext("Count",side=2,xpd=T,adj=0.5,padj=-3.7,las=0)
for (i in 2:ncol(pairend.meanqual)) {
  points(seq(0.5,nrow(pairend.meanqual)-0.5,1),log(pairend.meanqual[,i]+1),type="l",col=palette()[i]) 
}
# Mean base quality per position
plot(NULL,bty="n",xlim=c(0,nrow(pairend.meanposqual)),ylim=c(floor(min(pairend.meanposqual[,-1],na.rm=T)/10)*10,40),xlab="nucleotide position",ylab="Phred score")
for (i in 2:ncol(pairend.meanposqual)) {
  points(pairend.meanposqual[,1],pairend.meanposqual[,i],type="l",col=palette()[i]) 
}
# Length
lmax<-max(sum.pairend.length[,2])
posx<-unique(round(pairend.length[,1]/50))*50
posy<-ceiling(lmax/10^(nchar(lmax)-1))*10^(1:(nchar(lmax)-1))
barplot(log(sum.pairend.length[,2]+1),axes=F,border=NA,space=0,ylim=c(0,log(posy[length(posy)]+1)),col="limegreen",xlab="Read length")
axis(2,at=c(0,log(posy+1)),labels=c(0,sub("\\+0","\\+",format(posy,scientific=T))))
axis(1,at=posx-pairend.length[1,1],labels=posx)
mtext("Count",side=2,xpd=T,adj=0.5,padj=-3.7,las=0)
for (i in 2:ncol(pairend.length)) {
 points(seq(0.5,nrow(pairend.length)-0.5,1),log(pairend.length[,i]+1),type="l",col=palette()[i]) 
}
# Overlap
lmax<-max(sum.pairend.overlap[,2])
posx<-unique(round(pairend.overlap[,1]/25))*25
posy<-ceiling(lmax/10^(nchar(lmax)-1))*10^(1:(nchar(lmax)-1))
barplot(log(sum.pairend.overlap[,2]+1),axes=F,border=NA,space=0,ylim=c(0,log(posy[length(posy)]+1)),col="limegreen",xlab="Overlap region length")
axis(2,at=c(0,log(posy+1)),labels=c(0,sub("\\+0","\\+",format(posy,scientific=T))))
axis(1,at=posx,labels=posx)
mtext("Count",side=2,xpd=T,adj=0.5,padj=-3.7,las=0)
for (i in 2:ncol(pairend.overlap)) {
  points(seq(0.5,nrow(pairend.overlap)-0.5,1),log(pairend.overlap[,i]+1),type="l",col=palette()[i]) 
}
dev.off()


## Legend ##
pdf(paste(subp,"legend.pdf",sep="."),paper="a4",width=0,height=0,title=paste(subp,"sample legend"))
layout(matrix(c(1,2),ncol=1),heights=c(1,20))
par(mar=c(0,2,0,0),cex=1)
plot.new()
text(0.5,0.85,"Libraries legend",font=2,xpd=T,cex=1.4)
par(mar=c(0,0,0,0),cex=1)
plot.new()
legend("top",colnames(fwd.meanqual)[-1],lty=1,col=palette()[2:ncol(fwd.meanqual)],ncol=2,bty="n",lwd=4,y.intersp=1.5,text.width=c(0.4,0.4))
dev.off()
