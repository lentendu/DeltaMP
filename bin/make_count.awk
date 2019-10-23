BEGIN{
	l=split(S,h)
	printf "Representative_Sequence\ttotal"
	for (i=1;i<=l;i++){
		printf "\t%s",h[i]
	}
	printf "\n"
}
{
	if(NR==1){
		for (i in h){
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
			for (i in h){
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
	for (i in h){
		printf "\t%s", samp[h[i]]
	}
	printf "\n"
}
