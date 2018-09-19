library(plyr)
subp<-commandArgs()[7]
lib<-commandArgs()[8]
bin<-commandArgs()[9]
source(file.path(bin,"palette.R"))

# Sequence length in each sample
tmp_l<-scan(paste0(lib,".length"),what="",sep="\n")
len<-strsplit(tmp_l,"\t")
names(len) <- sapply(len, `[[`, 1)
len <- lapply(len, `[`, -1)

# Average sequence quality in each sample
tmp_q<-scan(paste0(lib,".meanqual"),what="",sep="\n")
qual<-strsplit(tmp_q,"\t")
names(qual) <- sapply(qual, `[[`, 1)
qual <- lapply(qual, `[`, -1)

# Average position quality in each sample
tmp_p<-scan(paste0(lib,".meanposqual"),what="",sep="\n")
pos<-strsplit(tmp_p,"\t")
names(pos) <- sapply(pos, `[[`, 1)
pos <- lapply(pos, `[`, -1)

# Plots
# For all samples
pdf(paste(subp,lib,"raw_reads_statistics.pdf",sep="."),paper="a4",width=0,height=0,title=paste(subp,lib,"raw reads statistics"))
layout(matrix(c(1:4),ncol=1),heights=c(1,4,4,4))
## Title
par(mar=c(0,1,0,0),cex=1)
plot.new()
text(0.5,0.85,"Raw reads statistics",font=2,xpd=T,cex=1.4)
text(0,0.3,paste0("Project: ",subp,"\nLibrary: ",lib),adj=c(0,NA),font=2,xpd=T,cex=1)
## length distribution
par(mar=c(4,4,3,1),mgp=c(2.5,0.8,0),cex=0.9)
hist(as.numeric(unlist(len)),bty="n",xlab="Length",ylab="n",breaks="Scott",col="brown2",border="brown2",main="Read length distribution")
## quality distribution
hist(as.numeric(unlist(qual)),bty="n",xlab="Phred score",ylab="n",breaks="Scott",col="steelblue2",border="steelblue2",main="Quality distribution")
## Quality over the sequence
mm<-max(laply(pos,length))
pos_tab<-as.matrix(ldply(pos,function(x) as.numeric(c(x,rep(0,mm-length(x)))))[,-1])
barplot(colSums(pos_tab)/colSums((pos_tab>0)*1),space=0,border=NA,xaxt="n",col="limegreen",ylim=c(0,40),ylab="Phred score",xlab="Sequence position",main="Quality distribution over the sequences")
axis(1,seq(0,1000,200),lwd=0.5,lwd.ticks=1)
## Nucleotide content
dev.off()

if (length(len)>1) {
	# For each samples
	pdf(paste(subp,lib,"sample_based_raw_reads_statistics.pdf",sep="."),paper="a4",width=0,height=0,title=paste(subp,lib,"sample based raw reads statistics"))
	layout(matrix(c(1,1,2,5,3,5,4,5),ncol=2,byrow=T),widths=c(4,1),heights=c(1,4,4,4))
	## Title
	par(mar=c(0,1,0,0),cex=1)
	plot.new()
	text(0.5,0.85,"Sample based raw reads statistics",font=2,xpd=T,cex=1.4)
	text(0,0.3,paste0("Project: ",subp,"\nLibrary: ",lib),adj=c(0,NA),font=2,xpd=T,cex=1)
	## length distribution
	par(mar=c(4,4,3,1),mgp=c(2.5,0.8,0),cex=0.9)
	color<-palette()[1:length(len)]
	length_tab<-llply(len,function(x) table(as.numeric(x)))
	mm<-laply(len,function(x) max(as.numeric(x)))
	plot(NULL,xlim=c(0,max(mm)),ylim=c(0,max(unlist(length_tab))),bty="n",xlab="Length",ylab="n",main="Read length distribution")
	for(i in 1:length(len)) {
	  points(as.numeric(names(length_tab[[i]])),length_tab[[i]],type="l",col=color[i])
	}
	## quality distribution
	qual_tab<-llply(qual,function(x) table(round(as.numeric(x))))
	plot(NULL,xlim=c(0,40),ylim=c(0,max(unlist(qual_tab))),bty="n",xlab="Phred score",ylab="n",main="Quality distribution")
	for(i in 1:length(qual)) {
	  points(as.numeric(names(qual_tab[[i]])),qual_tab[[i]],type="l",col=color[i])
	}
	## Mean quality per position for each sample
	mm<-laply(pos,function(x) length(as.vector(x)))
	plot(NULL,xlim=c(0,max(mm)),ylim=c(0,40),bty="n",xlab="Sequence position",ylab="Mean quality score",main="Quality distribution over the sequences")
	for(i in 1:length(pos)) {
	  points(seq(1,mm[i]),as.numeric(pos[[i]]),type="l",col=color[i])
	}
	## Legend
	par(mar=c(0,0,0,0),cex=1)
	plot.new()
	legend("center",names(pos),pch=16,col=color,bty="n")
	dev.off()
}

# Additionnal
## length vs. quality
#plot(laply(len,function(x) mean(as.numeric(x))),laply(qual,function(x) mean(as.numeric(x))),cex=1.5,pch=16,col=color,bty="n",xlab="Mean length (nt)",ylab="Mean quality (Phred score)")
