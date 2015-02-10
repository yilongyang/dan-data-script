#read subgenome assigned Dan data files list
#output one file for each LG containing the number of each category
open(IN,"C:\\Yang\\dan-sub-list.txt");
while(defined($line=<IN>)){
	chomp$line;
	$in=$line;
	$l=0;
	$in=~/-LG(\d)(.*?).txt-sub/;
	$match=$2;
	$match=~s/\s+/_/g;
	print "$match\n";
	if(!exists $file{$match}){
		push(@files,$match);
		$file{$match}++;
	}	
		
	open(DA,"$in")||die;
	while(defined($daline=<DA>)){
		chomp$daline;
		@temp=split("\t",$daline);
		if($temp[0]=~/EX/){
		#	next;
		}	
		$tag="$temp[1]";
		$tag=~s/\?/N/;
		$hlg=$temp[7];
		$hlg=~/LG(\d)/;
		$lg="LG".$1;
		$pos=$temp[6];
		if(!exists $alltag{$tag}){
			$alltag{$tag}++;
			push(@tags,$tag);
		}
		# if(!exists $hlgs{$hlg}){
			# push(@hmelgs,$hlg);
			# $hlgs{$hlg}++;
		# }
		$countkey="$tag-$lg-$match";
		$count{$countkey}++;
		#print "$in\t$lg\t$pos\t$countkey\t$count{$countkey}\n";
	}
	close DA;
}
#@matched=sort{lc($a) cmp lc($b)} @files;
for($i=1;$i<8;$i++){
	$lg="LG".$i;
	$out="Dan-$lg-all-count.txt";
	open(OUT,">>$out");
	
	foreach $tag(@tags){
		$header="";
		$number="";
		foreach $match(@files){
			$header=$header."$match\t";
			$countkey="$tag-$lg-$match";
			if(!exists $count{$countkey}){
				$count{$countkey}=0;
			}	
			$number=$number."$count{$countkey}\t";
		}
		print OUT "$tag\t$number\n";
	}
	print OUT "$lg\t$header\n";
	close OUT;
}
	
		
		
		
		