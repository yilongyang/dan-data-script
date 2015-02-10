#read octoploid Fiinumae Fvesca haplotype files
#mark msnps that come from one subgenome
#input sorted extended SNP file containing Fv, Fi, Octoploid haplotypes
open(DA,"C:\\Yang\\extend-haplotypes-octo-2-sorted.txt");
while(defined($line=<DA>)){
	chomp$line;
	@temp=split("\t",$line);
	$lg=$temp[0];
	$cpos=$temp[1];
	$mpos=$temp[3];
	$fvc=$temp[4];
	$fvc=~s/-/=/;
	$fii=$temp[5];
	$fii=~/^(.*?)-(\w)$/;
	$fiic=$1;
	$fiim=$2;
	$fiic=~s/-/=/;
	%fvsub=();
	%fiisub=();
	%funsub=();
	%fvsub2=();
	%fiisub2=();
	%funsub2=();
	$octohap="";
	if($fii=~/[A-Z]/ and $temp[6]=~/[A-Z]/ and $fiic!~$fvc){
		$substart=$cpos-1;
		$tag="N\tN\tN";
		$tag2="N\tN\tN";
		for($i=6;$i<scalar@temp;$i++){
			#print "octo haps\t$temp[$i]\n";
			$temp[$i]=~/\|(.*?)-(\w)/;
			$csnp=$1;
			$msnp=$2;
			#print "Check\t$temp[$i]\t$fvc\t$csnp\t$msnp\t$fvsub{$csnp}\n";
			if($csnp=~/$fvc/ and $fvc=~/$csnp/){
				if(!exists $fvsub{$csnp}){
					$fvsub{$csnp}=$msnp;
				}elsif(exists $fvsub{$csnp} and $fvsub{$csnp}!~$msnp){
					$tag=~s/^N\t/Y\t/;		# Fvesca is defined as having the same cSNP nucleotide, but any different mSNP nucleotide
				}
			}elsif($csnp=~/$fiic/ and $fiic=~/$csnp/){
				if(!exists $fiisub{$csnp}){
					$fiisub{$csnp}=$msnp;
					#print "first\t$csnp\t$msnp\t$fiisub{$csnp}\n";
				}elsif(exists $fiisub{$csnp} and $fiisub{$csnp}!~$msnp){
					#print "fii\t$csnp\t$msnp\t$fiisub{$csnp}\n";
					$tag=~s/\tN\t/\tY\t/;	# Fiinumae is defined as having the same cSNP nucleotide with Fiinumae sequence, but with any different mSNP nucleotide
				}							# cSNP must be different from Fvesca,if cSNP is an indel, it may be marked as both Fvesca and F.iinumae
			}else{
				if(!exists $funsub{$csnp}){
					$funsub{$csnp}=$msnp;
				}elsif(exists $funsub{$csnp} and $funsub{$csnp}!~$msnp){
					$tag=~s/\tN$/\tY/;	#Unknown is defined as having a cSNP different from both Fvesca and Fiinumae nucleotides,but with any different mSNP nucleotides.
				}
			}
			
			# if($csnp=~/$fiic/ and $fiic=~/$csnp/){
				# if(!exists $fiisub2{$csnp}){
					# $fiisub2{$csnp}=$msnp;
					# #print "first\t$csnp\t$msnp\t$fiisub{$csnp}\n";
				# }elsif(exists $fiisub2{$csnp} and $fiisub2{$csnp}!~$msnp){
					# #print "fii\t$csnp\t$msnp\t$fiisub{$csnp}\n";
					# $tag2=~s/\tN\t/\tY\t/;	# Fiinumae is defined as having the same cSNP nucleotide with Fiinumae sequence, but with any different mSNP nucleotide
				# }
			# }elsif($csnp=~/$fvc/ and $fvc=~/$csnp/){
				# if(!exists $fvsub2{$csnp}){
					# $fvsub2{$csnp}=$msnp;
				# }elsif(exists $fvsub2{$csnp} and $fvsub2{$csnp}!~$msnp){
					# $tag2=~s/^N\t/Y\t/;		# Fvesca is defined as having the same cSNP nucleotide, but any different mSNP nucleotide
				# }							# cSNP must be different from Fvesca,if cSNP is an indel, it may be marked as both Fvesca and F.iinumae
			# }else{
				# if(!exists $funsub2{$csnp}){
					# $funsub2{$csnp}=$msnp;
				# }elsif(exists $funsub2{$csnp} and $funsub2{$csnp}!~$msnp){
					# $tag2=~s/\tN$/\tY/;	#Unknown is defined as having a cSNP different from both Fvesca and Fiinumae nucleotides,but with any different mSNP nucleotides.
				# }
			# }
			
			
			
			
			
		}
		
		print "$tag\t$line\n";
		$count{$tag}++;
	}	
	
	
}
foreach $key(keys %count){
	print "$key\t$count{$key}\n";
}
	
			




				