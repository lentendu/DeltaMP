{
	if(NR==1) {
		id=$1
		sim=$2
		if(NF<4) {
			flag=0
		} else {
			flag=1
			ref[NR]=$3
			nf=split($4,line,",")
			maxnf=nf
			for(i=1;i<=nf;i++){
				grid[i,NR]=line[i]
			}
		}
	} else {
		ref[NR]=$3
		nf=split($4,line,",")
		nf>maxnf?maxnf=nf:maxnf
		for(i=1;i<=nf;i++){
			grid[i,NR]=line[i]
		}
	}
} END {
	printf "%s %s ",id,sim
	if(flag==0) {
		printf "%s", "NA"
		for(i=2;i<=8;i++) {
			printf ";%s", "NA"
		}
		printf "; %s\n","no_match"
	} else {
		if(FNR==1) {
			printf "%s(100)", grid[1,1]
			for(i=2;i<=nf;i++) {
				printf ";%s(100)" , grid[i,1]
			}
			printf "; %s\n",ref[1]
		} else {
			endref=ref[1]
			for(i=1;i<=maxnf;i++) {
				for(j=1;j<=FNR;j++) {
					if(length(grid[i,j])==0) {
						taxo[j]="NA"
					} else {
						taxo[j]=grid[i,j]
					}
				}
				rank[i]="NA"
				count=1
				prev=taxo[1]
				refall=ref[1]
				for (j=2;j<=FNR;j++){
					if(taxo[j]==prev) {
						count+=1
						if(ref[j-1]!=ref[j]) {
							refall=refall","ref[j]
						}
					} else {
						if(count/FNR*100>=cons) {
							rank[i]=prev
							stat[i]=count/FNR*100
							break
						} else {
							prev=taxo[j]
							count=1
							refall=ref[j]
						}
					}
				}
				if(count/FNR*100>=cons && rank[i]=="NA")
					rank[i]=prev
					stat[i]=count/FNR*100
				if(rank[i]!="NA")
					endref=refall
			}
			if(rank[i]=="NA") {
				printf "%s" , rank[1]
			} else {
				printf "%s(%.0f)" , rank[1], stat[1]
			}
			for(i=2;i<=nf;i++) {
				if(rank[i]=="NA") {
					printf ";%s" , rank[i]
				} else {
					printf ";%s(%.0f)" , rank[i], stat[i]
				}
			}
			printf "; %s\n", endref
		}
	}
}
