
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

%PE_type2col=("DEL","#c34343","DUP","#0375b6","INV","#a968b5","TRA","#f8a748","BND","#f8a748");

BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use My::Seq qw(GC_seq);
use My::VCF2; 


$readLen=150;


$samples=shift(@ARGV);
$control_BW_stem=shift(@ARGV);
$inRef=shift(@ARGV);
$projectDir=shift(@ARGV);
$inVCF=shift(@ARGV);

@Asamples=split(",",$samples);
print STDERR join(",",@Asamples)."\n";

$verbose=0;

print STDERR "# in vcf: $inVCF\n";


print STDERR "parse VCF...\n";
############################################################
############################################################
############################################################

my $VCFObj=My::VCF2->new($inVCF);

$VCFObj->add_header("##FORMAT=<ID=PAFDRA,Number=1,Type=Float,Description=\"# Percent of control samples with this event. A control sample has this event, if at least 90% of the event length is represented.   \">");


print ${$VCFObj->print_header};


while (my $V=$VCFObj->next_line){
  	
  	($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{"INFO"}{"END"});
	next if $cChr eq "hs37d5";
	next if $cChr eq "NC_007605";
		
	if ($$V{"INFO"}{"CNV"} == 0 or $$V{"INFO"}{"CHR2"} ne $$V{"CHROM"}){
		print ${$V->print_line}; next; 
	}
	
	if (!exists($$V{"INFO"}{"GC"})){
		$cReg="$cChr:$cSt-$cEn";
		$cSeq=""; open(GCIN,"samtools faidx $inRef $cReg |"); while(<GCIN>){next if $.==1; chomp; $cSeq.=$_}close(GCIN);
		$$V{INFO}{GC}=sprintf("%.0f",GC_seq(\$cSeq));
	}
	
  	print $_ if $verbose;
  	
  		 
	# GT:PE:SR:DRF:DBP:MQBP:FT:CNDR:CNP	
	
	foreach $cSample (@Asamples){
		
		########### calculate: Percent of control samples with this read depth ratio plus minus 0.1
		if(exists($$V{GTS}{DRA}) and exists($$V{GTS}{DRA}{$cSample}) ){
			
			$$V{GTS}{PAFDRA}{$cSample}=calcPercContr($cChr,$cSt,$cEn,$$V{GTS}{DRA}{$cSample});
		}
	}
	
	print ${$V->print_line};

# 	last if $c++ > 100;
}

sub calcPercContr{
	my $v=0;
	my ($cChr,$cSt,$cEn,$cDRA)=@_;
	
	($cSt,$cEn)=($cEn,$cSt) if $cSt>$cEn;
	
	$cLen=($cEn-$cSt+1);
	for($cExp=4; $cExp>=0; $cExp-- ){ # get right key for summarized data
		last if $cLen/(10**($cExp+3)) > 50;
	} $cExp=0 if $cExp<0;
	$cExpV=(10**($cExp+3));
	print "cLen: $cLen, cExp: $cExp, cExpV:$cExpV \n" if $v;	
	
	$caSt=int((($cSt-2)/$cExpV))+1; # skip partial interval, start from full 1k interval downstream, skip < 1-1001,  1001 -> 1, 1002-2001->2 (corresponding to start at position 2001)
	$caEn=int((($cEn-1)/$cExpV))-1; # same for end interval
	
	return (-1) if ($caSt>$caEn) or $caEn<0 or $caSt<0;
	return (-1) if !$cDRA or $cDRA<0;
	
	
	print "0: calcPercContr $cChr,$cSt,$cEn,$cDRA  $caSt:$caEn\n" if $v;
	###### if new chr load array
	if($visChr{$cChr}++ == 0){
		undef @covArr;
		if(! -e "$control_BW_stem/array/$cChr.tab" ){ 
			print STDERR "skipping chr: \"$cChr\" since $control_BW_stem/array/$cChr.tab does not exist\n";
		}else{
			print "read $cChr...\n" if $v;
			open(INC,"<$control_BW_stem/array/$cChr.tab") || die "does not exists $control_BW_stem/array/$cChr.tab";
			while(<INC>){
				chomp; @_=split("\t",$_);
				$cPos=($_[1]<1001)? 0:int((($_[1]-1)/1000));
				$covArr[0][$cPos]=[@_[2..$#_]];
	# 				print STDERR "$_[1] $cPos  $_  \n".join("\t",@$covArr[$cPos])."\n";exit if $ccc++>100;
			}close(INC);
			
			for my $cL (1..4){ # summerize the data in intervals of 10k, 100k, 1M, 10M
				$lL=$cL-1;
				for(my $i=0; ($i+9)<=$#{$covArr[$lL]}; $i+=10 ){
					undef @tmpArr;
					for my $j ($i..($i+9)){ # for each coordinate
						for(my $k=0; $k<=$#{$covArr[$lL][$j]}; $k++ ){ # for each value
							$tmpArr[$k]+=$covArr[$lL][$j][$k];
						}
# 						print "j: $j ".int($i/10)." ".join("\t", @{$covArr[$lL][$j]} )."\n" if $v;
					}
					for(my $k=0; $k<=$#{$covArr[$lL][$j]}; $k++ ){ # compute the average
						$tmpArr[$k]=round($tmpArr[$k]/10,1);
					}
# 					print "cL: $cL ".int($i/10)." ".join("\t",@tmpArr)."\n" if $v;
					$covArr[$cL][int($i/10)]=[@tmpArr];					
				}
			}
			
		}
	}	
	
	return (-1) if ! $covArr[$cExp];
	
	print "  cDRA:$cDRA \n"  if $v;
	undef @sample2segmentCountAll;
	undef @sample2segmentCountEvent;
		
	###### calc average per colum
	undef @tmpArr;
	for $cPos ($caSt..$caEn){ # for each pos
		for $cSmpl (0..$#{$covArr[$cExp][$cPos]}){
			
			push @{$tmpArr[$cSmpl]}, $covArr[$cExp][$cPos][$cSmpl];
			
		}
	}
	
	for $cSmpl (0..$#tmpArr){
		undef @tmpArr2;
		
		print "raw\tsample$cSmpl\t".join("\t",@{$tmpArr[$cSmpl]})."\n"  if $v;	
		for $cPos (0..$#{$tmpArr[$cSmpl]}){ # if possible smooth by averaging with neighboring 2 or 1
	
			my $cVal=0;my $cC=0;
			if( ($cPos-2)>=0 and ($cPos+2)<=$#{$tmpArr[$cSmpl]} ){  foreach (($cPos-2)..($cPos+2)){ $cV=$tmpArr[$cSmpl][$_]; next if $cV<0; $cC++; $cVal+=$cV;   }   $cVal= $cC>0 ? round($cVal/$cC,1):(-1);    }
			elsif( ($cPos-1)>=0 and ($cPos+1)<=$#{$tmpArr[$cSmpl]} ){  foreach (($cPos-1)..($cPos+1)){ $cV=$tmpArr[$cSmpl][$_]; next if $cV<0; $cC++; $cVal+=$cV;   }   $cVal= $cC>0 ? round($cVal/$cC,1):(-1);    }
			else{ $cVal=$tmpArr[$cSmpl][$cPos]; }
			$tmpArr2[$cPos]=$cVal;
							
		}
		print "smoothed\tsample$cSmpl\t".join("\t",@tmpArr2)."\n"  if $v;	
		
		for $cPos (0..$#tmpArr2){ # remove salt
	
			if( $tmpArr2[$cPos]<=1.2 and $tmpArr2[$cPos]>=0.8 ){
				my $cVal=0;my $cC=0;my $cCAll=0;  
				foreach (($cPos-6)..($cPos+6)){ next if $_<0 or $_>$#tmpArr2; $cV=$tmpArr2[$_]; next if $cV<0; $cCAll++; next if $cV>=0.8 and $cV<=1.2; $cC++; $cVal+=$cV }
				if($cC/$cCAll>=0.8){  $tmpArr2[$cPos]= $cC>0 ? round($cVal/$cC,1):(-1);  }			
			}
		}
		
		print "desalted\tsample$cSmpl\t".join("\t",@tmpArr2)."\n"  if $v;	
		
		for $cPos (0..$#tmpArr2){ # remove salt	
			$cVal=$tmpArr2[$cPos];
			next if $cVal<0; 
			if($cDRA>=1.2 and $cVal>=1.2 ){ $sample2segmentCountEvent[$cSmpl]++; }
			elsif($cDRA<=0.8 and $cVal<=0.8 ){ $sample2segmentCountEvent[$cSmpl]++; }
			elsif($cDRA>=0.8 and $cDRA<=1.2  and  $cVal>=0.8 and $cVal<=1.2 ){ $sample2segmentCountEvent[$cSmpl]++; }
			$sample2segmentCountAll[$cSmpl]++;					
		}
			
	}
		
	$samples90pcSame=0;
	$samplesAll=0;
	
	for $cSmpl (0..$#sample2segmentCountAll){
		$sample2segmentCountEvent[$cSmpl]+=0;
		$sample2segmentCountAll[$cSmpl]+=0;
		
		print "1: cSmpl:$cSmpl event:$sample2segmentCountEvent[$cSmpl] all:$sample2segmentCountAll[$cSmpl]\n" if $v;
		next if $sample2segmentCountAll[$cSmpl]==0;
		$samplesAll++;
		next if $sample2segmentCountEvent[$cSmpl]==0;
		next if $sample2segmentCountEvent[$cSmpl]/$sample2segmentCountAll[$cSmpl]<0.9;
		$samples90pcSame++;
	}	
	
	print "2: cFreq:$cFreq samplesAll:$samplesAll samples90pcSame:$samples90pcSame\n" if $v;
	
	return -1 if $samplesAll==0;
	$cFreq=$samples90pcSame/$samplesAll;
	$RL=($cFreq==0)? 0:int(log10($cFreq)*(-1))+2;
	$cFreq=round($cFreq,$RL); 	
	#print STDERR "".join(",",@covArrAvg)."\n".join(",",@subCovArrAvg)."\n$cChr,$cSt,$cEn $caSt,$caEn $lowLim,$upLim:".$$V{GTS}{DRA}{$cSample}." ".scalar(@subCovArrAvg)."/".scalar(@covArrAvg)." $cFreq\n";
	#print STDERR "$cChr,$cSt,$cEn $cLen $caSt,$caEn $lowLim,$upLim:".$cDRA." ".scalar(@subCovArrAvg)."/".scalar(@covArrAvg)." $cFreq\n";
	return $cFreq;
	
}



	


sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}
sub log10 { my $n = shift;  return log($n)/log(10); }
























