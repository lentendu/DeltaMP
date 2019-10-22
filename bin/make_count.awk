BEGIN{
	gsub(";","\t",S)
	l=split(S,h)
	print "Representative_Sequence\ttotal\t"S
}
{
	printf "%s\t%s",$1,$2
	for(i=1;i<=l;i++){
		if($3~h[i]":"){
			printf "\t%s", gensub(".*"h[i]":([0-9]*).*","\\1",$3)
		} else {
			printf "\t%s",0
		}
	}
printf "\n"
}
