#assign final subgenome to each mSNP
open(DA,"C:\\Yang\\extend-haplotypes-octo-2-sorted-subgenome.txt");
while(defined($line=<DA>)){
	chomp$line;
	@temp=split("\t",$line);
	$fv=$temp[0];
	$fii=$temp[1];
	$fun=$temp[2];
	$lg=$temp[3];
	$msnp=$temp[6];
	$msnpkey="$lg-$msnp";
	$tag="$fv-$fii-$fun";
	$countkey="$msnpkey-$tag";
	#print "$line\n";
	
	if($msnpkey!~/$lastmsnpkey/){
		$l=0;
		foreach $key(sort{$count{$b}<=>$count{$a}} keys %count){
			$num=$count{$key};
			$key=~/LG(\d)-(\d+)-(.*?)$/;
			#print "Key\t$key\n";
			$subtag=$3;
			$tags=$tags."$subtag\t$num\t";
			$l++;
			
		}
		%count=();
		$lastmsnpkey=~s/-/\t/g;
		# 	$tags=~s/-/\t/g;
		print "$l\t$lastmsnpkey\t$tags\n";
		$tags="";
		
	}
	$count{$countkey}++;
	#print "$countkey\t$count{$countkey}\n";
	$lastmsnpkey=$msnpkey;
}
	
	