BEGIN{
	l=split(S,h,"#")
}
{
	if(NR==1){
		for (i=1;i<=l;i++){
			samp[h[i]]=0
		}
		prev=$1
		sum=$4
		samp[$3]=$4
	} else {
		if($1==prev){
			sum+=$4
			samp[$3]=$4
		} else {
			printf "%s\t%s", prev,sum
			for (i=1;i<=l;i++){
				printf "\t%s", samp[h[i]]
				samp[h[i]]=0
			}
			printf "\n"
			prev=$1
			sum=$4
			samp[$3]=$4
		}
	}
} END {
	printf "%s\t%s", prev,sum
	for (i=1;i<=l;i++){
		printf "\t%s", samp[h[i]]
	}
	printf "\n"
}
