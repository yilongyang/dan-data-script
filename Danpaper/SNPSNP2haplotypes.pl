#read IStraw90 SNP-SNP list
#At the corresponding positions, find haplotypes from F. iinumae and octoploid sam files
@samfiles=(
"C:\\Yang\\Rosbreed\\Fragaria_iinumae-HD-2004-15_vs_fvesca_v1.1_pseudo.fna.all.nodups.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa_DovxCam_F2_34_CRAG_vs_fvesca_v1.1_pseudo.fna.all.nodups.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Holiday-71006_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Korona_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Winter_Dawn_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Fanella_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa_vs_fvesca_v1.1_psuedo.fna.all.sorted.nodups.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Emily_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
"C:\\Yang\\Rosbreed\\Fragaria_x_ananassa-Sweet_Charlie_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
);

open(IN,"C:\\Yang\\Rosbreed.7764snpsnp.11112012.csv");
while(defined($line=<IN>)){
	chomp$line;
	@linetemp=split(",",$line);
	$center=$linetemp[1];
	@snptemp=split(qw(\|),$center);
	$lg=$snptemp[1];
	$csnppos=$snptemp[2];
	$msnppos=$snptemp[5];
	$csnpkey="$lg-$csnppos";
	$msnpkey="$lg-$msnppos";
	$csnp{$csnpkey}=$center;
	$msnp{$msnpkey}=$center;
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
		$cigar=~s/(\d+)([A-Z])/$1$2-/g;
		@cigars=split("-",$cigar);
		$substart=0;
		$lgstart=$start-1;
		$nextsnp=$snp[$i];
		@nexttemp=split("-",$nextsnp);
		$nextlg=$nexttemp[0];
		$nextpos=$nexttemp[1];
		$distance=$nextpos-$start;
		# if($lg!~/LG1/ or $start>100000){
			# last;
		# }	
		#print "$lg\t$start\t$nextlg\t$nextpos\n";
		#print "$distance\t$lg\t$nextlg\t$i\tNext snp\n";
		if($distance>-1 and $distance < 100 and $lg=~/$nextlg/  and $mapq>20){
		
			#print "\n$seq\n"
			#print "File $l\t$lg\t$start\t$nextlg\t$nextpos\n";
			foreach $cigars(@cigars){
				if($cigars=~/I|S/){
					$cigars=~/(\d+)/;
					$d=$1;
					$substart=$substart+$d;
				}elsif($cigars=~/D/){
					$cigars=~/(\d+)/;
					$d=$1;
					#$del=substr($lgseq{$lg}, $lgstart,$d);
					#print "$lg, $lgstart,$d, DEL\t$del\n";
					#$readseq=$readseq.$del;
					$lgstart=$lgstart+$d;
				}elsif($cigars=~/M/){
					$cigars=~/(\d+)/;
					$d=$1;
					
					$match=substr($seq, $substart,$d);
					$qmatch=substr($qual,$substart,$d);
					#$readseq=$readseq.$match;
					@matchtemp=split("",$match);
					@qmatchtemp=split("",$qmatch);
					$haplotype="csnp-msnp";
					#print "$readname\t$substart\t$d\n";
					$lgpos=$lgstart;
					for($k=0;$k<scalar@matchtemp;$k++){
						$baseq=$basequality{$qmatchtemp[$k]};
						#print "$readname\t$seq\n$qual\n$matchtemp[$k]\t$k\t$baseq\n";
							
						$lgpos++;
						$basekey="$lg-$lgpos";
						#print "$basekey\n";
						if($baseq<30){
							next;
						}
						if(exists $csnp{$basekey}){
							$haplotype=~s/csnp/$matchtemp[$k]/;
							#print "CSNP\t$base\n";
						}elsif(exists $msnp{$basekey}){
							$haplotype=~s/msnp/$matchtemp[$k]/;
							$msnpkey=$basekey;
							#print "MSNP\t$base\n";
						}
						
						
						if($haplotype!~/snp/){
							$haplotypekey="$msnpkey-$haplotype";
							#print "Haplotype\t$haplotype\n";
							if(!exists $hap{$haplotypekey}){
								$hap{$haplotypekey}++;
								if($l==1){
									$fiihap{$haplotypekey}++;
									$haplotypes{$msnpkey}=$haplotypes{$msnpkey}."=$l|$haplotype";
								#}elsif($l==2){
								#	$haplotypes{$msnpkey}=$haplotypes{$msnpkey}."=$l|$haplotype";
								}else{
									$octohap{$haplotypekey}++;
									$haplotypes{$msnpkey}=$haplotypes{$msnpkey}."\t$l|$haplotype";
								}			
								print "$l\t$msnpkey\t$haplotypes{$msnpkey}\n";
							}else{	
								$hap{$haplotypekey}++;
								$fiihap{$haplotypekey}++;
								$octohap{$haplotypekey}++;
								#print "$haplotypekey\t$octohap{$haplotypekey}\n";
								
							}	
							last;				
							
						}
					}	
					$substart=$substart+$d;
					$lgstart=$lgstart+$d;
					#print "OUT\n";
				}
				#print   "OUT Cigar\n";
			}
		}elsif($distance<0 && $lg=~/$nextlg/  ||  $distance>0 && $lg!~/$nextlg/){
			$i++;
			#print "$distance\t$lg\t$nextlg\t$i\tNext snp\n";
		}elsif($lg!~/$nextlg/){
			next;
		}
		
		#print "Next read\n";	
	}
	if($l==1){
		%hap=();	
	}	
	
	close SAM;
}
$outfile="snpsnp-haplotypes-5.txt";
open(OUT,">>$outfile");	
foreach $snp(@snp){
	@temp=split("\t",$haplotypes{$snp});
	$fiihaplotype=shift@temp;
	@fiitemp=split(qw(=1\|),$fiihaplotype);
	shift@fiitemp;
	foreach $fiitype (@fiitemp){
		$haplotypekey="$snp-$fiitype";
		if($fiihap{$haplotypekey}>1){
			$fiihaps=$fiihaps."=1|$fiitype";	
				
			}
		
		#$haplotypes{$snp}=~s/=1\|$fiitype/=1\|$fiitype<$fiihap{$haplotypekey}/;
	}	
	foreach $octotype(@temp){
		@octosplit=split(qw(\|),$octotype);
		$num=shift@octosplit;
		$foctotype=$octosplit[0];
		$haplotypekey="$snp-$foctotype";
		#print "-$octotype\t-$foctotype\n";
		
		#print "Octo count\t$haplotypekey\t$octohap{$haplotypekey}\n";
		if($octohap{$haplotypekey}>1){
			
			$foctohaps=$foctohaps."$num|$foctotype\t";
		}	
		
	}	
		
	print OUT "$snp\t$msnp{$snp}-iinumae$fiihaps\t$foctohaps\n";	
	$fiihaps="";
	$foctohaps="";
}	
