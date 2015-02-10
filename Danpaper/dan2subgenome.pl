#find markers used in Dan data
#read assiged markers from all extended list
open(IN,"C:\\Yang\\extend-haplotypes-ocgto-2-sorted-subgenome-final.txt");
while(defined($line=<IN>)){
	chomp$line;
	@temp=split("\t",$line);
	$extendkey="$temp[1]-$temp[2]";
	if($temp[0]==1){
		$extend{$extendkey}=$temp[3];
	}	
}
close IN;
#read assignmed markers from snpsnp list	
open(IN,"C:\\Yang\\snpsnp-haplotypes-5-subgenome.txt");
while(defined($line=<IN>)){
	chomp$line;
	@temp=split("\t",$line);
	$snpkey=$temp[3];
	$snp{$snpkey}="$temp[0]-$temp[1]-$temp[2]";
}
close IN;
open(IN,"C:\\Yang\\dan-list.txt");
$out="C:\\Yang\\dan-stat.txt";
open(OTH,">>$out");
while(defined($line=<IN>)){
	chomp$line;
	$out=$line."-sub.txt";
	open(OUT,">>$out");
	open(DA,"$line")||die;
	$agree=0;
	$disagree=0;
	$snpsnp=0;
	$extend=0;
	while(defined($daline=<DA>)){
		chomp$daline;
		@temp=split("\t",$daline);
		$snpkey="$temp[2]-$temp[3]";
		$snpsub=$snp{$snpkey};
		$snpsub=~s/\?/N/g;
		if(exists $snp{$snpkey}){
			if($snpsub=~/$extend{$snpkey}/ and $extend{$snpkey}=~/\w/){
				print OUT "A\t$snp{$snpkey}\t$extend{$snpkey}\t$daline\n";
				#print "A\t$snp{$snpkey}\t$extend{$snpkey}\n";
				$agree++;
				if(!exists $asnp{$snpkey}){
					$asnp{$snpkey}++;
				}	
			}elsif($snpsub!~/$extend{$snpkey}/ and $extend{$snpkey}=~/\w/){
				print OUT "B\t$snp{$snpkey}\t$extend{$snpkey}\t$daline\n";
				#print "B\t$snp{$snpkey}\t$extend{$snpkey}\n";
				$disagree++;
				if(!exists $dsnp{$snpkey}){
					$dsnp{$snpkey}++;
				}	
			}else{
				print OUT "C\t$snp{$snpkey}\t$extend{$snpkey}\t$daline\n";
			}	
			$snpsnp++;
			$adsnp{$snpkey}++;
		}elsif(exists $extend{$snpkey}){
			print OUT "EX\t$extend{$snpkey}\t\t$daline\n";
			$extend++;
			if(!exists $exsnp{$snpkey}){
				$exsnp{$snpkey}++;
			}	
		}
		
	}
	print OTH "$snpsnp\t$agree\t$disagree\t$extend\t$line\n";
	close DA;
	close OUT;
}
	$asnpsum=scalar keys%asnp;
	print OTH "Agree markers\t$asnpsum\n";
	$dsnpsum=scalar keys%dsnp;
	print OTH "Disagree markers\t$dsnpsum\n";
	$adsnpsum=scalar keys%adsnp;
	print OTH "All SNPSNP markers\t$adsnpsum\n";
	$exsnpsum=scalar keys%exsnp;
	print OTH "Extended markers\t$exsnpsum\n";