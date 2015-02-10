#read snpsnp-Da filtered snps and dertermine subgenome
open("IN","C:\\Yang\\snpsnp-haplotypes-5.txt");
$out1="C:\\Yang\\snpsnp-haplotypes-5-subgenome.txt";
$out2="C:\\Yang\\snpsnp-haplotypes-5-other.txt";
open(OUT,">>$out1");
open(OTH,">>$out2");
while(defined($line=<IN>)){
	chomp$line;
	$subtag="N\tN\tN";
	@temp=split("\t",$line);
	$center=$temp[1];
	$center=~/mSNP\|LG(\d)\|(\d+)\|(\w)\/(\w)\|LG(\d)\|(\d+)\|(\w)\/(\w)-iinumae(.*?)$/;
	#print "$2-$3-$5-$6\n";
	$fvcsnp=$3;
	$altcsnp=$4;
	$refm=$7;
	$altm=$8;
	$vesca1="$fvcsnp-$refm";
	$vesca2="$fvcsnp-$altm";
	$un1="$altcsnp-$refm";
	$un2="$altcsnp-$altm";
	
	$iinumae=$9;
	#print "Fvesca\t$vesca1\t$vesca2\nUnknown\t$un1\t$un2\n";
	shift@temp;
	shift@temp;
	$octohap=join("",@temp);
	
	if($iinumae=~/[A-Z]/){
	
		if($octohap=~/$vesca1/ and $octohap=~/$vesca2/){
			$subtag=~s/^N/Y/;
		}	
		#print "Fiinumae\t$iinumae\n";
		@fitemp=split(qw(=1\|),$iinumae);
		shift@fitemp;
		$fvfi=0;
		foreach $fi(@fitemp){
			print "Fi\t$fi\n";
			$fi=~/(\w)-(\w)/;
			$iinumaesnp=$1;
			##print "$iinumaesnp\n";
			$fi1="$iinumaesnp-$refm";
			$fi2="$iinumaesnp-$altm";
			if($iinumaesnp=~/$fvcsnp/){
				$fvfi++;
			}	
					
			if($octohap=~/$fi1/ and $octohap=~/$fi2/){
				#if($iinumaesnp=~/$fvcsnp/){
					$subtag=~s/\tN\t/\tY\t/;
				if($iinumaesnp=~/$altcsnp/ and $subtag=~/^Y/){
					$subtag="cSNP-conflict";
					print OTH "$subtag\t$line\n";
					last;
				}elsif($iinumaesnp!~/$fvcsnp/ and $iinumaesnp!~/$altcsnp/ and $subtag=~/^Y/){
					$subtag="other-conflict";
					print OTH "$subtag\t$line\n"; #$octohap\n$fi1\t$fi2\n";
					last;
				}	
			}
		}
		if($octohap=~/$un1/ and $octohap=~/$un2/ and $fvfi>0 ){
			
			if($subtag=~/^Y/){
				$subtag="cSNP-conflict";
				print OTH "$subtag\t$line\n";
			}else{
				$subtag=~s/\tN$/\tY/;
			}	
		}	
		$subtag=~s/Y\tN\tN/Y\tN\t?/;
		$subtag=~s/N\tY\tN/N\tY\t?/;
		$subtag=~s/Y\tN\tY/Y\tN\t?/;
		if($subtag!~/conflict/){
			print OUT "$subtag\t$line\n";
		}	
	}else{
		%check=();
		%badcc=();
		for($i=0;$i<scalar@temp;$i++){
			$temp[$i]=~/(\w)-(\w)/;
			$cc=$1;
			$mm=$2;
			#print OTH "$temp[$i]\n";
			if(!exists $check{$cc}){
				$check{$cc}=$mm;
			}elsif($check{cc}!~/$mm/){
				$badcc{$cc}++;
			}	
		}
		$num=scalar keys%badcc;
		if($num>1){
			print OTH "other-conflict\t$line\n";
		}else{
			print OTH "$subtag\t$line\n";
		}
		
	}
	

}
	
				
				
				
				
				
				
				
		