
# Author: Andre Minoche, a.minoche@garvan.org.au

use Getopt::Long;

GetOptions ("s=s" => \$stepS,
			"r=s"  => \$region,
		    "o=s"  => \$outStem,
		    "f=s"  => \$fastaIn,
		    "b=s"  => \$bamIn,
		    "a"  => \$append
)  or die("Error in command line arguments\n");


# ($stepS,$region,$outStem,$fastaIn,$bamIn)=@ARGV;


die "can not parse region: $region " if $region !~ /^(.+):([0-9,]+)-([0-9,]+)$/;

($chr,$st,$en)=($1, $2, $3) if ($region =~ /^(.+):([0-9,]+)-([0-9,]+)$/);

print STDERR "region: $chr,$st,$en\n";

$verbose=1;

if ($append){ 
	open(OUT,">>$outStem.wig");
}else{
	open(OUT,">$outStem.$chr.wig");
}

open(FASTA,"samtools faidx $fastaIn $chr |") || die "can not open samtools faidx $ARGV[2] $chr |";
open(SMQ, "samtools mpileup --reference $fastaIn -s -r $chr $bamIn | cut -f 1,2,7 |") || die "can not open samtools mpileup -r $chr $ARGV[3] |";

print OUT join(" ",("fixedStep","chrom=$chr","start=$st","step=$stepS","span=$stepS"))."\n";


$pG=0;
while($f=<FASTA>){ chomp($f); next if $.==1; $fBuf.=$f; $fBufLen=length($fBuf); $fBufPos=0;

	while(($fBufPos+$stepS)<=$fBufLen){
	
	 $sumMQ=0;$nonN_count=0;
# 	 $verbStr1="",$verbStr2="" if $verbose;
	 for $p ($fBufPos..($fBufPos+$stepS-1)){
	  $pG++;
	  if(substr($fBuf,$p,1) ne "N"){$nonN_count++; $sumMQ+=getMQ($pG);}
# 	  if($verbose){ $verbStr1.=substr($fBuf,$p,1);   $verbStr2.=getMQ($pG).", "; }
	 }
# 	 print "verbStr\tlenBuf:".length($fBuf).">= $stepS, $verbStr1 $verbStr2 $fBuf\n" if $verbose;
	 if($stepS==1){$cAvgDepth=($nonN_count>0)? $sumMQ:(-1);}
	 else{$cAvgDepth=($nonN_count>=($stepS/2))? sprintf("%.0f",($sumMQ/$nonN_count)):(-1);}
# 	 print join("\t",($pG,$nonN_count,$sumMQ,$cAvgDepth))."\n" if $verbose;
# 	 print OUT $cAvgDepth." $pG $fBufPos ".($fBufPos+$stepS-1)." ".$fBufLen."\n"; 
	 print OUT $cAvgDepth."\n"; 
	 $fBufPos+=$stepS;
	}
	$fBuf=substr($fBuf,$fBufPos);
# 	last if $pG>=100000;
} 


close(FASTA); 
close(SMQ);


sub getMQ{ my ($fP)=@_; #1       10001   N       14      ^IT^!T^!T^!T^!T^(T^2T^JT^!T^!T^=T^!T^!T^!t      FAFAAAAAFAK<A<  I!!!!(2J!!=!!!

 if(!$dL){$dL=<SMQ>; chomp($dL); ($pC,$pP,$pMQ)=split("\t",$dL); $avgMQ=avgMQ(\$pMQ); print join("\t",("mpileup1",$dL))."\n" if 0; }
 
 while($pP<$fP){ return 0 unless $dL=<SMQ>; chomp($dL); ($pC,$pP,$pMQ)=split("\t",$dL); $avgMQ=avgMQ(\$pMQ);  print join("\t",("mpileup2",$dL))."\n" if 0; }
 
 if($pP>$fP){return (-2);} # for regions without mapping quality return -2 instead of 0 for depth
 elsif($pP==$fP){return $avgMQ;}
 else{die "s! $dL"}
}


sub avgMQ{
	my ($MQs)=@_;
	my $MQ_count=length($$MQs);
	my $MQ_sum=0;
	my $MQ_avg=(-2);

	return $MQ_avg if length($$MQs)==0;

	for ($i=0; $i<=($MQ_count-1);$i++){
		$MQ_sum+=(ord(substr($$MQs,$i,1))-33);
	}

	$MQ_avg=sprintf("%.0f",$MQ_sum/$MQ_count);

	return $MQ_avg;
}











