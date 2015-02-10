#read snp list from DaxMo data
#At the corresponding positions, find haplotypes from F. iinumae sam files
@samfiles=(
"Rosbreed/Fragaria_iinumae-HD-2004-15_vs_fvesca_v1.1_pseudo.fna.all.nodups.bam.sam",

);
open(IN,"marker-list-sorted.txt")||die;
while(defined($line=<IN>)){
	chomp$line;
	@linetemp=split("\t",$line);
	$lg=$linetemp[0];
	$msnppos=$linetemp[1];
	$msnpkey="$lg-$msnppos";
	$msnp{$msnpkey}=$linetemp[2];
	push(@snp,$msnpkey);
	#print "$msnpkey\n";
}
close IN;	
$i=0;
$l=0;
$quality=qw(!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~);

@quality=split("",$quality);
for($kk=0;$kk<=scalar@quality;$kk++){
	$basequality{$quality[$kk]}=$kk;
}
foreach $file(@samfiles){
	open(SAM,"$file")||die;
	$l++;
	$i=0;
	while(defined($samline=<SAM>)){
		@samlinetemp=split("\t",$samline);
		$readname=$samlinetemp[0];
		$lg=$samlinetemp[2];
		$start=$samlinetemp[3];
		$mapq=$samlinetemp[4];
		$cigar=$samlinetemp[5];
		$seq=$samlinetemp[9];
		$qual=$samlinetemp[10];
		#print "$lg-$start\n";
		$samline=~/MD:Z:(.*?)(\s+|\n)/;
		$md=$1;
		$cigar=~s/(\d+)([A-Z])/$1$2-/g;
		@cigars=split("-",$cigar);
		$substart=0;
		$lgstart=$start-1;
		$lgmdstart=$start;
		$nextsnp=$snp[$i];
		@nexttemp=split("-",$nextsnp);
		$nextlg=$nexttemp[0];
		$nextpos=$nexttemp[1];
		$distance=$nextpos-$start;
		$m=0;
		@readexsnp=();
		#print "$lg\t$start\t$nextlg\t$nextpos\n";
		#print "$distance\t$lg\t$start\t$nextlg\t$nextpos\n";
		if($distance>-1 and $distance < 100 and $lg=~/$nextlg/  and $mapq>20){
			#print "MD\t$md\n";
			$md=~s/(\d+)([A-Z])/$1-$2-/g; #extract reference allele from the MD string
			@md=split("-",$md);
			foreach $mdx(@md){
				if($mdx=~/\d+/){
					$lgmdstart=$lgmdstart+$mdx;
				}elsif($mdx=~/\^(\w+)/){
					$del=$1;
					$len=length$del;
					$lgmdstart=$lgmdstart+$len;
				}elsif($mdx=~/\w+/){
					$exsnpkey="$lg-$lgmdstart";
					if(!exists $exsnp{$exsnpkey}){
						$exsnp{$exsnpkey}=$mdx;
						#print "exSNP\t$exsnpkey\n";
					}	
					$lgmdstart=$lgmdstart+length$mdx;
				}	
			}
			#print "\n$seq\n"
			#print "File $l\t$lg\t$start\t$nextlg\t$nextpos\n";
			for($k=0;$k<scalar@cigars;$k++){
				if($cigars[$k]=~/I|S/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;
					$ins=substr($seq,$substart,$d);
					#$inskey="$lg-$lgstart-I-$ins";
					#push(@readexsnp,$inskey); #
					#$depth{$inskey}++;	# this hash count the frequency of this type of insertion
					$exsnpkey="$lg-$lgstart";
					#if(!exists $exsnp{$exsnpkey}){
					#	$exsnp{$exsnpkey}="I";
						
					#}	
								
					$substart=$substart+$d;
				}elsif($cigars[$k]=~/D/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;
					#$delkey="$lg-$lgstart-D-$d";
					#push(@readexsnp,$delkey);
					#$depth{$delkey}++; #count
					#$exsnpkey="$lg-$lgstart";
					#if(!exists $exsnp{$exsnpkey}){
					#	$exsnp{$exsnpkey}="D";
					#}
					
					$lgstart=$lgstart+$d;					
				}elsif($cigars[$k]=~/M/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;					
					$match=substr($seq, $substart,$d); #extract the aligned sequence fragment
					$qmatch=substr($qual,$substart,$d); #extract the corresponding quality score
					#$readseq=$readseq.$match;
					@matchtemp=split("",$match);
					@qmatchtemp=split("",$qmatch);
					#print "$readname\t$substart\t$d\n";
					$lgpos=$lgstart;
					for($mk=0;$mk<scalar@matchtemp;$mk++){
						$baseq=$basequality{$qmatchtemp[$mk]};
					
						$lgpos++;
						$basekey="$lg-$lgpos";
						#print "SNP\t$basekey\n";
						if($baseq<30){
							next;
						}
						
						if(exists $msnp{$basekey}){
							#print "MSNP\t$base\n";
							$msnpbase=$matchtemp[$mk];
							$msnppos=$basekey;
							$m=1;
						}elsif(exists $exsnp{$basekey}){
							$readbase="$lg-$lgpos-$matchtemp[$mk]";
							push(@readexsnp,$readbase);
							$depth{$readbase}++; #cout;
							#print "not mSNP\t$readbase\n";
						}	
					}
					$substart=$substart+$d;
					$lgstart=$lgstart+$d;
					
				}
				
			}	
			if($m>0){	
				foreach $excsnp(@readexsnp){
					
					$excsnp=~/LG(\d)-(\d+)/;
					$fvCSNP="LG".$1."-$2";
					$snppair="$exsnp{$fvCSNP}-$excsnp=$msnpbase";
					if(!exists $check{$snppair}){
						$haplotype{$msnppos}=$haplotype{$msnppos}."|$snppair";
						$check{$snppair}=$snppair;
						print "Found\t$msnppos\t$snppair\n";
					}
				}	
			}				
		}elsif($distance<0 && $lg=~/$nextlg/  ||  $distance>0 && $lg!~/$nextlg/){
			$i++;
			#print "$distance\t$lg\t$nextlg\t$i\tNext snp\n";
		}elsif($lg!~/$nextlg/){
			next;
		}		
		#print "Next read\n";	
		
		
	
	}
	close SAM;
}
$outfile="extend-haplotypes-iinumae-1.txt";
open(OUT,">>$outfile");	
foreach $snp(@snp){
	@temp=split(qw(\|),$haplotype{$snp});
	foreach $haplotype(@temp){
		$haplotype=~/LG(.*?)=/;
		$exsnp="LG".$1;
		if($depth{$exsnp}>1){
			$fiihap=$fiihap."|$haplotype";
		}
	}		
	print OUT "$snp\t"."iinumae$fiihap\n";	
	$fiihap="";
}	
