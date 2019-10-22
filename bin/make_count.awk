BEGIN{
	gsub(" ","\t",S)
	l=split(S,h)
	print "Representative_Sequence\ttotal\t"S
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
