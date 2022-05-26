
# Author: Andre Minoche, a.minoche@garvan.org.au



use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";


use Bio::DB::Big;



$bwType=shift(@ARGV);
$project=shift(@ARGV);
$cChr=shift(@ARGV);
$outWigStem=shift(@ARGV);

$digits=0;
$digits=1 if $bwType eq "mq";

@bigwigs=@ARGV;

@Asamples=split(",",$samples);



print STDERR "open the bw files...\n";
############################################################
our %bwObj;

foreach $cQBWpre (@bigwigs){ 
	foreach $cT (("mq","q0","q20")){
			
		$cQBW=$cQBWpre.".$cT.bw";
		die " $cQBW ! exists" if (! -e $cQBW);
		$bwObj{$cQBWpre}{$cT} = Bio::DB::Big->open($cQBW);
	}	
	print STDERR "# sample: $cQBWpre, bw found\n";
}


print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
############################################################

@meanCov;

$c=0;
foreach $cQBWpre (@bigwigs){
	$sumCov=0;
	for(1..22){
		$sumCov+=getCov($_,20000001,30000000,"q0",58,$cQBWpre); #($qChr,$qSt,$qEn,$covQ,$qMQ,$cSample)=@_;
	}
	$meanCov[$c]=sprintf("%.2f",($sumCov)/22);
	print STDERR "average coverage $cSample is: ".$meanCov[$c]."\n";
	$c++;
}

# print STDERR join(",",@meanCov)."\n";

print STDERR "open wig files...\n";
############################################################
$c=0;
%fh;
foreach $cQBWpre (@bigwigs){
	
	$cQBW=$cQBWpre.".$bwType.bw";
	print STDERR " bigWigToWig -chrom=$cChr $cQBW /dev/stdout \n";
	open($fh{$c},"bigWigToWig -chrom=$cChr $cQBW /dev/stdout | ") || die " bigWigToWig -chrom $cChr $cQBW /dev/stdout \n";
	$c++;	
}




print STDERR "create average wig...\n";
############################################################


$NrSamples=@bigwigs;
$NrSamplesMid=int($NrSamples/2);
@v;

open(OUT,">$outWigStem/control10.$bwType.$cChr.wig") || die " con not write to \n ";
if($bwType eq "q0"){ open(DEV,">$outWigStem/control10.dv.$cChr.wig") || die " con not write \n "; }
$continue=1;
while($continue==1){

	for $i (0..($NrSamples-1)){
		
		$cFH=$fh{$i};
		if($v[$i]=<$cFH>){ chomp($v[$i]); }else{ $continue=0 }
		
		
		
		if(  substr($v[$i],0,1) eq "f"  ){ # if header print and get next number
		
			print OUT $v[$i]."\n" if $i==0; # print fixedStep chrom=1 start=1 step=1 span=1
			print DEV $v[$i]."\n" if($bwType eq "q0") and $i==0;
 			if($v[$i]=<$cFH>){ chomp($v[$i]); }else{ $continue=0 }
			
		}
# 		print  "$v[$i] $meanCov[$i]\n" if $v[$i]>0;
		$v[$i]=($v[$i]<0)? (-1):$v[$i]/$meanCov[$i]; # normalize
# 		print  "$i: $v[$i]\n";# if $v[$i]>0;
		
	}
	last if $continue==0;
	
	@v=sort {$a <=> $b} @v;
	
	if($NrSamples%2==0){		
		$cOutVal=sprintf("%.".$digits."f",(($v[($NrSamplesMid-1)]+$v[$NrSamplesMid])/2)*100);
	}else{
		$cOutVal=sprintf("%.".$digits."f",$v[$NrSamplesMid]*100);	
	}
	
	$cOutVal=(-1) if $cOutVal<0;
	print OUT $cOutVal."\n";
	print DEV aad(\@v)."\n" if($bwType eq "q0");
	
}
close(OUT);
close(DEV) if($bwType eq "q0");

# 0 1 2 3 4 5 6 7 8 9
# 1 2 3 4 5 6 7 8 9

sub getCov{
	my ($qChr,$qSt,$qEn,$covQ,$qMQ,$cSample)=@_;
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	
	$qSt=0 if $qSt<0;
	$qEn=$chr2len{$qChr} if $qEn > $chr2len{$qChr};
	return (-1) if $qSt>=$qEn;
	
	my $qLen=($qEn-$qSt+1);
	my ($sumCov,$fBins)=(0,0);
    my $qBins=100;

    $qBins=$qLen if ($qLen/$qBins)<1;

	my $aCov_stat = $bwObj{$cSample}{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	my $aMQ_stat = $bwObj{$cSample}{mq}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean') if $qMQ>0;	
	
	for ($i=0; $i<=$#$aCov_stat;$i++){		
		
      $cMQ=$qMQ>0 ? ${$aMQ_stat}[$i]:0;
      $cCov=${$aCov_stat}[$i];

	  next if $qMQ>0 and $cMQ<$qMQ;
	  next if $cCov<0;	  
	  $sumCov+=$cCov; $fBins++;
	  		
	}

	
	if($fBins==0 or $fBins/$qBins<0.33333){ return (-1) }
	else{ return sprintf("%.2f",$sumCov/$fBins) }
}


sub aad {

($cA)=@_;

 my $count=0;
 my $sum=0;
 foreach (@$cA) { next if $_== (-1); $count++;   $sum+=$_  }
 return (-1) if $count==0;
 my $avg=$sum/$count;

 my $sumDev=0;
 foreach (@$cA) {  next if $_== (-1); $sumDev+=abs($_-$avg);  }

return sprintf("%.2f",$sumDev/$count);

}


















