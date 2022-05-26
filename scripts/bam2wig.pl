
# Author: Andre Minoche, a.minoche@garvan.org.au

use Getopt::Long;

GetOptions ("s=s" => \$stepS,
			"q=s"  => \$MapQCutOff,
		    "r=s"  => \$region,
		    "o=s"  => \$outStem,
		    "f=s"  => \$fastaIn,
		    "b=s"  => \$bamIn,
		    "a"  => \$append
)  or die("Error in command line arguments\n");


# ($stepS,$MapQCutOff,$region,$outStem,$fastaIn,$bamIn)=@ARGV;


$MapQCutOff=~ s/q//g;

die "can not interprete MQ $MapQCutOff " if $MapQCutOff !~ /^[0-9]+$/;


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
open(SDEPTH, "samtools depth --reference $fastaIn -d 500000 -Q $MapQCutOff -r $chr $bamIn |") || die "can not open samtools depth -r $chr $ARGV[3] |";

print OUT join(" ",("fixedStep","chrom=$chr","start=$st","step=$stepS","span=$stepS"))."\n";


$pG=0;
while($f=<FASTA>){ chomp($f); next if $.==1; $fBuf.=$f; $fBufLen=length($fBuf); $fBufPos=0;

	while(($fBufPos+$stepS)<=$fBufLen){
	
	 $sumDepth=0;$nonN_count=0;
# 	 $verbStr1="",$verbStr2="" if $verbose;
	 for $p ($fBufPos..($fBufPos+$stepS-1)){
	  $pG++;
	  if(substr($fBuf,$p,1) ne "N"){$nonN_count++; $sumDepth+=getDepth($pG);}
# 	  if($verbose){ $verbStr1.=substr($fBuf,$p,1);   $verbStr2.=getDepth($pG).", "; }
	 }
# 	 print "verbStr\tlenBuf:".length($fBuf).">= $stepS, $verbStr1 $verbStr2 $fBuf\n" if $verbose;
	 if($stepS==1){$cAvgDepth=($nonN_count>0)? $sumDepth:(-1);}
	 else{$cAvgDepth=($nonN_count>=($stepS/2))? sprintf("%.1f",($sumDepth/$nonN_count)):(-1);}
# 	 print join("\t",($pG,$nonN_count,$sumDepth,$cAvgDepth))."\n" if $verbose;
# 	 print OUT $cAvgDepth." $pG $fBufPos ".($fBufPos+$stepS-1)." ".$fBufLen."\n"; 
	 print OUT $cAvgDepth."\n"; 
	 $fBufPos+=$stepS;
	}
	$fBuf=substr($fBuf,$fBufPos);
# 	last if $pG>=100000;
} 


close(FASTA); 
close(SDEPTH);

sub getDepth{ my ($fP)=@_;
 if(!$dL){$dL=<SDEPTH>; chomp($dL); ($dC,$dP,$dD)=split("\t",$dL); print join("\t",("depth1",$dL))."\n" if 0; }
 
 while($dP<$fP){ return 0 unless $dL=<SDEPTH>; chomp($dL); ($dC,$dP,$dD)=split("\t",$dL);  print join("\t",("depth2",$dL))."\n" if 0; }
 
 if($dP>$fP){return 0;}
 elsif($dP==$fP){return $dD;}
 else{die "s! $dL"}
}

