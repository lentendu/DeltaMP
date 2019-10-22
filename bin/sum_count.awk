{
	if(NR==1){
		seq=$1
		prev=$2
		sum=$3
		for(i=4;i<=NF;i++){
			a[i]=$i
		}
	} else {
		if($2==prev){
			sum+=$3
			for(i=4;i<=NF;i++){
				a[i]+=$i
			}
		} else {
			printf "%s\t%s",seq,sum
			for(i=4;i<=NF;i++){
				printf "\t%s",a[i]
				a[i]=$i
			}
			printf "\n"
			seq=$1
			prev=$2
			sum=$3
		}
	}
} END {
	printf "%s\t%s",seq,sum
	for(i=4;i<=NF;i++){
		printf "\t%s",a[i]
	}
	printf "\n"
}
