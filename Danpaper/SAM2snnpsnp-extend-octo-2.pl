#read extend Fiinumae SNP list file
#for each row extract one or more cSNP-mSNP pairs
#search for haplotypes of those cSNP-mSNP pairs in a set of stawberry cultivars



@samfiles=(
"Rosbreed/Fragaria_x_ananassa_DovxCam_F2_34_CRAG_vs_fvesca_v1.1_pseudo.fna.all.nodups.realigned.bam.sam",

"Rosbreed/Fragaria_x_ananassa-Holiday-71006_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa-Korona_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa-Winter_Dawn_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa-Fanella_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa_vs_fvesca_v1.1_psuedo.fna.all.sorted.nodups.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa-Emily_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
"Rosbreed/Fragaria_x_ananassa-Sweet_Charlie_vs_fvesca_v1.1_pseudo.fna.all.sorted.nodups.fixed.realigned.bam.sam",
);
open(IN,"extend-haplotypes-iinumae-1.txt");
while(defined($line=<IN>)){
	chomp$line;
	$line=~/^LG(\d)-(\d+)\tii/;
	$lg="LG".$1;
	$mpos=$2;
	$msnp="$lg-$mpos";
	
	push(@snp,$msnp);
	@linetemp=split("\t",$line);
	$msnp{$msnp}=$linetemp[1];
	@ctemp=split(qw(\|),$linetemp[1]);
	shift@ctemp;
		
	foreach $csnpinfor(@ctemp){
		if($csnpinfor!~/\w/){
			last;
		}	
		$csnpinfor=~/(\w)-LG(\d)-(\d+)-(.*?)=(\w)/;
		$csnpkey="LG"."$2-$3";
		$csnpfv{$csnpkey}=$1;
		$fihap="$4-$5";
		$hapfikey="$csnpkey-$msnp";
		if(!exists $hapfi{$hapfikey}){
			$hapfi{$hapfikey}=$fihap;
		}	
		#print "Hapfikey\t$hapfikey\t$fihap\n";
		
	}	
	
	
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
		$m=0;
		@readexsnp=();
		#print "$lg\t$start\t$nextlg\t$nextpos\n";
		#print "$distance\t$lg\t$nextlg\t$i\tNext snp\n";
		if($distance>-100 and $distance < 100 and $lg=~/$nextlg/  and $mapq>20){
			#print "\n$seq\n"
			#print "File $l\t$lg\t$start\t$nextlg\t$nextpos\n";
			for($k=0;$k<scalar@cigars;$k++){
				if($cigars[$k]=~/I|S/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;
					# $ins=substr($seq,$substart,$d);
					# $inskey="$lg-$lgstart-I=$ins";
					# if($cigars[$k]=~/I/){
						
						# $exsnpkey="$lg-$lgstart";
						# if(exists $csnpfv{$exsnpkey}){
							# push(@readexsnp,$inskey);
							
						
						# }
					# }		
								
					$substart=$substart+$d;
				}elsif($cigars[$k]=~/D/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;
					# $delkey="$lg-$lgstart-D=$d";
					
					# $exsnpkey="$lg-$lgstart";
					# if(exists $csnpfv{$exsnpkey}){
						# push(@readexsnp,$delkey);
						
					# }
					
					$lgstart=$lgstart+$d;					
				}elsif($cigars[$k]=~/M/){
					$cigars[$k]=~/(\d+)/;
					$d=$1;					
					$match=substr($seq, $substart,$d);
					$qmatch=substr($qual,$substart,$d); 
					#$readseq=$readseq.$match;
					@matchtemp=split("",$match);
					@qmatchtemp=split("",$qmatch);
					#print "$readname\t$substart\t$d\n";
					$lgpos=$lgstart;
					for($mk=0;$mk<scalar@matchtemp;$mk++){
						$baseq=$basequality{$qmatchtemp[$mk]};
						$lgpos++;
						$basekey="$lg-$lgpos";
						if($baseq<30){
							next;
						}
						#print "SNP\t$basekey-$base\t$readname\n";
						if(exists $msnp{$basekey}){
							#print "MSNP\t$base\n";
							$msnpbase=$matchtemp[$mk];
							$msnppos=$basekey;
							$m++;
						}elsif(exists $csnpfv{$basekey}){
							$readbase="$lg-$lgpos-$matchtemp[$mk]";
							push(@readexsnp,$readbase);
							#print "SNP added\t$readbase\n";
						}	
					}
					$substart=$substart+$d;
					$lgstart=$lgstart+$d;
					
				}
				
			}	
			if($m>0){	
				foreach $readexsnp(@readexsnp){
					@temp=split("-",$readexsnp);
					$haplotypekey="$temp[0]-$temp[1]-$msnppos";
					$snppair="$readexsnp=$msnppos-$msnpbase";
					$count{$snppair}++;
					$readexsnp=~/-(\d+)-(.*?)$/;
					$exsnp=$2;
					$hap="$exsnp-$msnpbase";
					#print "Haplotypekey\t$haplotypekey\n";
					if(!exists $check{$snppair}){
						$haplotype{$haplotypekey}=$haplotype{$haplotypekey}."\t$l|$hap";
						$check{$snppair}++;
						print "Found\tMSNP$msnppos\t$snppair\n";
					}
				}	
			}				
		}elsif($distance<-100 && $lg=~/$nextlg/  ||  $distance>100 && $lg!~/$nextlg/){
			$i++;
			#print "$distance\t$lg\t$nextlg\t$i\tNext snp\n";
		}elsif($lg!~/$nextlg/){
			
			next;
			#print "Change CHR\n";
		}		
		#print "Next read\n";	
		
		
	
	}
	close SAM;
}
$outfile="extend-haplotypes-octo-2.txt";
open(OUT,">>$outfile");	
foreach $cmsnp(keys %haplotype){
	$cmsnp=~/^LG(\d)-(\d+)-(.*?)-(\d+)/;
	$csnp="LG"."$1-$2";
	$msnp="$3-$4";
	@temp=split("\t",$haplotype{$cmsnp});
	foreach $hap(@temp){
		$hap=~/(\d)\|(\w)-(\w)/;
		$num=$1;
		$csnpbase=$2;
		$msnpbase=$3;
		$snppair="$csnp-$csnpbase=$msnp-$msnpbase";
		if($count{$snppair}>1){
			$octohaplotype=$octohaplotype."$hap\t";
		}
	}	
	print OUT "$cmsnp\t$csnpfv{$csnp}\t$hapfi{$cmsnp}\t$octohaplotype\n";
	$octohaplotype="";
}
	
	
	
	
	
	
