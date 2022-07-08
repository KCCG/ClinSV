
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";


%PE_type2col=("DEL","#c34343","DUP","#0375b6","INV","#a968b5","TRA","#f8a748","BND","#f8a748");

%cat2type=(
"1p"=>"INS1", "2m"=>"INS1", "p"=>"INS1",
"2p"=>"INS2", "1m"=>"INS2", "m"=>"INS2",

"1pU:2pD"=>"DEL", "2mU:1mD"=>"DEL", "pU:mD"=>"DEL",
"1p:2p"=>"BND1", "2m:1m"=>"BND1", "p:m"=>"BND1",

"1mU:2mD"=>"DUP", "2pU:1pD"=>"DUP", "mU:pD"=>"DUP", "pD:mU"=>"DUP",
"1m:2m"=>"BND2", "2p:1p"=>"BND2", "m:p"=>"BND2",

"1pU:2mD"=>"INV1", "2mU:1pD"=>"INV1", "pU:pD"=>"INV1",
"1p:2m"=>"BND1", "2m:1p"=>"BND1", "p:p"=>"BND1",

"1mU:2pD"=>"INV2", "2pU:1mD"=>"INV2", "mU:mD"=>"INV2",
"1m:2p"=>"BND2", "2p:1m"=>"BND2", "m:m"=>"BND2");


use Bio::DB::Big;
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use My::Seq qw(GC_seq);
use File::Basename;

##### v4 
# - use bigwig bins, and 
# - add read depth ratio for MQ=60 portion, if portion > 30%
# - convert DEL or dups not passing the depth filter to BND, do not exclude

$readLen=150;
$inBED=shift(@ARGV);
$cSample=shift(@ARGV);
$inRef=shift(@ARGV);
$inRefStyle=shift(@ARGV);
$projectDir=shift(@ARGV);
$jobfs=shift(@ARGV);
$S_SV_control_numSamples=shift(@ARGV);
$S_control_bw_folder=shift(@ARGV);

if ($inRefStyle =~ /chr/){
	$chrPf='chr';
	$Y_chr="chrY";
	$X_chr="chrX";
	$M_chr="chrM";
}else{
	$chrPf='';
	$Y_chr="Y";
	$X_chr="X";
	$M_chr="MT";
}


%tabixC;
$tabixC{SR} = new Tabix(-data => "$jobfs/SR.brkp.gz",-index => "$jobfs/SR.brkp.gz.tbi");
$tabixC{PE} = new Tabix(-data => "$jobfs/PE.brkp.gz",-index => "$jobfs/PE.brkp.gz.tbi");

$verbose=1;

if($inBED=~/^(.+)\.(bed|txt)(.gz)*$/){ $outBED=$1 } $outBED.=".f1.bed";

print STDERR "# running filterCNVnator-v2.pl\n";
print STDERR "# in BED: $inBED\n";
print STDERR "# out BED: $outBED\n";


print STDERR "open splitters discordants as tabix...\n";
############################################################
$cIN=$projectDir."/SVs/$cSample/lumpy/bed/".$cSample.".splitters.brkp.gz"; die " $cIN ! exists" if (! -f $cIN);
`tabix -f -s 1 -b 2 -e 2  $cIN`;
$txSR = new Tabix(-data => "$cIN",-index => "$cIN.tbi");
open(IN1, " tabix --list-chroms $cIN | ") || die " tabix --list-chroms $cIN |  not found";
while(<IN1>){ chomp; $txSR_chroms{$_}++ }close(IN1);


$cIN=$projectDir."/SVs/$cSample/lumpy/bed/".$cSample.".discordants.brkp.gz"; die " $cIN ! exists" if (! -f $cIN);
`tabix -f -s 1 -b 2 -e 2  $cIN`;
$txPE = new Tabix(-data => "$cIN",-index => "$cIN.tbi");
open(IN1, " tabix --list-chroms $cIN | ") || die " tabix --list-chroms $cIN |  not found";
while(<IN1>){ chomp; $txPE_chroms{$_}++ }close(IN1);



print STDERR "open the bw files...\n";
############################################################
our %bwObj;


foreach $cT (("mq","q0","q20")){
			
	$cQBW=$projectDir."/alignments/$cSample/bw/".$cSample.".$cT.bw";
	print $cQBW."\n";
	die " $cQBW ! exists" if (! -f $cQBW);
	$bwObj{$cT} = Bio::DB::Big->open($cQBW);
}	
print STDERR "# sample: $cSample, bw found\n";

$bwObj{mq} = Bio::DB::Big->open($S_control_bw_folder."/aMQ.bw");
$bwObj{sd} = Bio::DB::Big->open($S_control_bw_folder."/popCovStdev.bw");


print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
############################################################


#read in chr2len
open(IN1, "<$inRef.fai") || die "can not open $inRef";
while(<IN1>){ chomp; @_=split("\t",$_); $chr2len{$_[0]}=$_[1]; }close(IN1);


print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
############################################################
%meanCov;
$sumCov=0;
for(1..22){
	$sumCov+=getCov("$chrPf$_",20000001,30000000,"q0",58,$cSample);
	print STDERR "$chrPf$_ $sumCov \n";
}

foreach $cChr ((keys %chr2len,"A")){ $meanCov{$cChr}=sprintf("%.2f",($sumCov)/22); }

# determine X and Y coverage and round to 0, 0.5, 1, 1.5, 2, ... times average coverage
$meanCov{$X_chr}=sprintf("%.0f",(getCov($X_chr,1,$chr2len{$X_chr},"q0",58,$cSample)/$meanCov{A})/0.5)*0.5*$meanCov{A};
$meanCov{$Y_chr}=sprintf("%.5f", (getCov($Y_chr,6641419,10079253,"q0",0,$cSample)+getCov($Y_chr,13800704,23668908,"q0",0,$cSample))/2 ); # two mq 60 region one of each Y arm 
print STDERR "Y:".$meanCov{$Y_chr}."\n";
$meanCov{$Y_chr}=0 if $meanCov{$Y_chr} == (-1); # if -1 then because it is not covered	
$meanCov{$Y_chr}=sprintf("%.0f",($meanCov{$Y_chr}/$meanCov{A})/0.5)*0.5*$meanCov{A};

# MT coverage
$meanCov{$M_chr}=sprintf("%.2f",(getCov($M_chr,1,$chr2len{$M_chr},"q0",58,$cSample)));	

print STDERR "average coverage $cSample autosome:".$meanCov{A}.", X:".$meanCov{$X_chr}.", Y:".$meanCov{$Y_chr}.", MT:".$meanCov{$M_chr}."\n";




# print STDERR "read chromosome length...\n";
# ############################################################
# %refID2Len;
# open (INS,"<$inRef.fai") || die "can not read $inRef.fai";
# while(<INS>){chomp; @_=split("\t",$_); $refID2Len{$_[0]}=$_[1];}
# close(INS);


print STDERR "parse BED...\n";
############################################################
############################################################
############################################################


$allVariantsinBED=0;
open(IN, "<$inBED") || die "$outBED nicht gefunden";
open(OUT, ">$outBED") || die "$outBED nicht gefunden";

while (<IN>){ #deletion        1:71801-73300   1500    0.449022        0.0282577       4.59695 1       1       0.394366        0.0278622

  	print $_ if $verbose;
	chomp; @_=split("\t",$_);
	($cChr,$cSt,$cEn)=split(/[\-:]/,$_[1]);
	($cType,$cFreq,$cPval,$cQ0F)=($_[0],$_[3],$_[4],$_[9]);
	$cLen=($en-$st+1);
	$cType=uc(substr($cType,0,3));
	next if $cChr eq "hs37d5";
	next if $cChr eq "NC_007605";
	
	print "chromosome: $cChr\n";
	print "Skipping: $cChr\n" and next if $cChr !~ /(^chr1[0-9]$)|(^chr2[0-2]$)|(^chr[1-9XYM]$)|(^1[0-9]$)|(^2[0-2]$)|^[1-9XY]$|(^MT$)/; 
	next if uc($cChr) eq $Y_chr and $meanCov{$Y_chr}==0 and $cType eq "DEL"; # ignore deletion on Y chromosome if Y has coverage == 0
	
	# split the CNVnator call at N regions if applicable
	undef @Segments;
	$pcValid=splitSeg($cChr,$cSt,$cEn,"q0",\@Segments);
	print "### next because >90% Ns\n" if $pcValid<0.1 and $verbose;
	next if $pcValid<0.1;
	
	# MQ of CNV
	undef %h;
	map { ($cSt,$cEn)=@$_; $h{"$cSt,$cEn"}{MQ}=getMQSD($cChr,$cSt,$cChr,$cEn,"mq");  } @Segments;
	map { ($cSt,$cEn)=@$_; $h{"$cSt,$cEn"}{PCSD}=getMQSD($cChr,$cSt,$cChr,$cEn,"sd");  } @Segments;
	
	# MQ of flanking 10% of CNV 
	map { ($cSt,$cEn)=@$_; $h{"$cSt,$cEn"}{MQBP}=calcMQBP($cChr,$cSt,$cChr,$cEn);  } @Segments;
	
	# depth of flanking 10% of CNV 
	map { ($cSt,$cEn)=@$_; $h{"$cSt,$cEn"}{DBP}=calcDBP($cChr,$cSt,$cEn);  } @Segments;
	
	# HQ read depth and % of CNV
	foreach  (@Segments){
		($cSt,$cEn)=@$_;
		($cDRF,$cDRA)=(-1,-1);
		calcDRF($cChr,$cSt,$cEn);
		$h{"$cSt,$cEn"}{DRF}=$cDRF;
		$h{"$cSt,$cEn"}{DRA}=$cDRA;
	}
	
	# get flanking SRs and PEs		
	map {
		($cSt,$cEn)=@$_;
		
		if( $h{"$cSt,$cEn"}{DRA}>50 ){
		
			foreach $cInfo (("SR","PE","CSR","CPE","PAFSU","SU")){
				$h{"$cSt,$cEn"}{$cInfo}=(-1);
			}
		}else{		
		
			$h{"$cSt,$cEn"}{SR}=concordantSRPEs($cChr,$cSt,$cEn,$cType,$txSR);  
			$h{"$cSt,$cEn"}{PE}=concordantSRPEs($cChr,$cSt,$cEn,$cType,$txPE); 
			$h{"$cSt,$cEn"}{CSR}=concordantSRPEs($cChr,$cSt,$cEn,$cType,$tabixC{SR}); 
			$h{"$cSt,$cEn"}{CPE}=concordantSRPEs($cChr,$cSt,$cEn,$cType,$tabixC{PE});
		
			$cSUC=($h{"$cSt,$cEn"}{CSR}+$h{"$cSt,$cEn"}{CPE});
			$cSU=$h{"$cSt,$cEn"}{SR}+$h{"$cSt,$cEn"}{PE};
			$h{"$cSt,$cEn"}{SU}=$cSU;
			$h{"$cSt,$cEn"}{PAFSU}=($cSU>0)? roundA(  ($cSUC/$cSU)/$S_SV_control_numSamples   ):(-1);  
		
			print "--> SR:".$h{"$cSt,$cEn"}{SR}.", PE:".$h{"$cSt,$cEn"}{PE}.", CSR:".$h{"$cSt,$cEn"}{CSR}.", CPE:".$h{"$cSt,$cEn"}{CPE}.", en:$cEn PAFSU:".$h{"$cSt,$cEn"}{PAFSU}."  \n" if $verbose;
		}
		
	} @Segments; 
		
	

	
	# print
	$cCol=$PE_type2col{$cType};
	
	foreach  (@Segments){
		($cSt,$cEn)=@$_;
		$cLen2=$cEn-$cSt+1;
		$cSeg="$cSt,$cEn";
		
		$cReg="$cChr:$cSt-$cEn";
		$cSeq=""; open(GCIN,"samtools faidx $inRef $cReg |"); while(<GCIN>){next if $.==1; chomp; $cSeq.=$_}close(GCIN);
		$h{$cSeg}{GC}=sprintf("%.0f",GC_seq(\$cSeq));
		
		# assign PASS
# 		if ($h{$cSeg}{MQBP}<20 or $h{$cSeg}{GC}>90 or $h{$cSeg}{GC}<10){
# 			$h{$cSeg}{FT}="FAIL";
# 		}else{
			$h{$cSeg}{FT}="PASS";
# 		}
		
		$CNVCount++;
		$gffTags="name=CNVnator_$CNVCount;TOOL=CNVnator;SAMPLE=$cSample;LOCATION=$cReg;SVTYPE=$cType;SVLEN=$cLen2;FREQ=$cFreq;PVAL=$cPval;Q0=$cQ0F;";
		foreach $cT (sort keys %{$h{$cSeg}}){ $gffTags.="$cT=".$h{$cSeg}{$cT}.";" if exists($h{$cSeg}{$cT}); }
		print OUT join("\t",($cChr,($cSt-1),$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen2,"0"))."\n";
		print join("\t",($cChr,($cSt-1),$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen2,"0"))."\n" if $verbose;
	
	}

}

close(IN);close(OUT);



sub splitSeg{ #splitSeg($cChr,$cSt,$cEn,"q0",\@Segments);
	my ($qChr,$qSt,$qEn,$covQ,$segsX)=@_;
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	$qSt=0 if $qSt<0;
	$qEn=$chr2len{$qChr} if $qEn > $chr2len{$qChr};
	return (-1) if $qSt>=$qEn;
	
	my $qLen=($qEn-$qSt+1);
	my ($sumCov,$fBins)=(0,0);
    my $qBins=100;
    $qBins=$qLen if ($qLen/$qBins)<1;
    
    $prevCov=(-1); $segC=(-1);
    
# 	undef @aCov;
# 	@aCov = $bwObj{$covQ}->features(-seq_id=>$qChr,-start=>$qSt,-end=>$qEn,-type=>"bin:$qBins");

# 	print "$qChr,$qSt,$qEn,$covQ,$cSample,$qLen, bins: $qBins \n";    
    
# 	for ($i=0; $i<=$#aCov;$i++){
# 	  my ($cSt,$cEn,$cCov)  = ($aCov[$i]->start,$aCov[$i]->end,$aCov[$i]->mean);
	  
	#print STDERR "splitSeg $qChr, $qSt, $qEn, $qBins\n";  
	my $aCov_stat = $bwObj{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	$bin_len=round($qLen/$qBins,0);
	
	for ($i=0; $i<=$#$aCov_stat;$i++){
	  
	  $cCov=${$aCov_stat}[$i];
	  $cSt=$qSt+$i*$bin_len;
	  $cEn=$cSt+$bin_len;
	  #print "       splitSeg $cSt $cEn \n";
	  
	  my $cLen = ($cEn-$cSt+1);

	  
# 	  print "c: $cSt,$cEn,$cLen: cCov:$cCov segC:$segC ";
	  if ($cEn >= $qEn){
		 if($prevCov>=0 and $segC>=0){ $$segsX[$segC][1]=$qEn }
		 last;
	  }
	 
	  if($cCov>=0){ # non N
		if($prevCov<0 or $i==0 ){ $segC++; $$segsX[$segC][0]=$cSt; } # start new segment
		if($segC>=0){ $$segsX[$segC][1]=$cEn; $$segsX[$segC][1]=$qEn if $i==$#aCov; }
		$fBins++;
		$sumCov+=$cCov;
	  }
	  $prevCov=$cCov;
# 	 print "     -> st,en: $$segsX[$segC][0]-$$segsX[$segC][1]\n";
	 
	}
# 	print " \n\n";
# 	map { print "seg $_, $$segsX[$_][0], $$segsX[$_][1] \n" } 0..$segC;
	
# 	print "f: $fBins ".sprintf("%.1f",$fBins/$qBins)."\n";
	return sprintf("%.1f",$fBins/$qBins);
}


sub getMQSD{
	my ($qChr,$qSt,$qChr2,$qEn,$MQorSD)=@_;
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	$qSt=0 if $qSt<0;
	$qEn=$chr2len{$qChr} if $qEn > $chr2len{$qChr};
	return (-1) if $qChr ne $qChr2 or $qSt>=$qEn;
	
	my $qLen=($qEn-$qSt+1);
	my ($sumMQ,$fBins)=(0,0);
    my $qBins=100;

    $qBins=$qLen if ($qLen/$qBins)<1;   
         
	undef @aMQ;
	#print STDERR  "getMQSD $qChr,$qSt,$qChr2,$qEn,$qLen, qMQ: $qMQ, bins: $qBins, MQorSD:$MQorSD\n";
# 	@aMQ = $bwObj{$MQorSD}->features(-seq_id=>$qChr,-start=>$qSt,-end=>$qEn,-type=>"bin:$qBins");	
#     
# 	for ($i=0; $i<=$#aMQ;$i++){		
# 	  my ($cSt,$cEn,$cMQ)  = ($aMQ[$i]->start,$aMQ[$i]->end,$aMQ[$i]->mean);
# 	  my $cLen = ($cEn-$cSt+1); 
# # 	  print "c: $cSt,$cEn,$cLen: $cCov,$cMQ ".($cMQ<$qMQ or $cCov<0)."\n"; 
# 	  next if $cMQ<0;
# 	  $sumMQ+=$cMQ; $fBins++;
# 	}
	
	my $aMQ_stat = $bwObj{$MQorSD}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	foreach my $cMQ (@{$aMQ_stat}) {
		#printf("%f\n", $s);
		next if $cMQ<0;
		$sumMQ+=$cMQ; $fBins++;
	}
	
# 	print "f: $fBins ".sprintf("%.1f",$sumMQ/($fBins+0.0000001))."\n";
	
	if($fBins==0 or $fBins/$qBins<0.33333){ return (-1) }
	else{ $rN=($MQorSD eq "mq")? 1:2; return round($sumMQ/$fBins,$rN) }
	
}

sub calcMQBP {

	($cChr,$cSt,$cChr2,$cEn)=@_;
	$cLen=($cEn-$cSt+1);
	$fl=int($cLen*0.1);
	
	$cMQBP1=getMQSD($cChr,($cSt-$fl),$cChr,($cSt-1),"mq");
	$cMQBP2=getMQSD($cChr,($cEn+1),$cChr,($cEn+$fl),"mq");
	$cMQBP=sprintf("%.0f",avgIfGr0($cMQBP1,$cMQBP2));	
	return $cMQBP;	
}

sub calcDBP {

	($cChr,$cSt,$cEn)=@_;
	$cLen=($cEn-$cSt+1);
	$fl=int($cLen*0.1);
	
	$cDBP1=getCov($cChr,($cSt-$fl),($cSt-1),"q0",0,$cSample);
	$cDBP2=getCov($cChr,($cEn+1),($cEn+$fl),"q0",0,$cSample);
	$cDBP=sprintf("%.0f",avgIfGr0($cDBP1,$cDBP2));	
	return $cDBP;	
}

sub calcDRF {

	($cChr,$cSt,$cEn)=@_;
	$cLen=($cEn-$cSt+1);
	
	# for flanking only rely on high MQ, else take global mean
	$cDFL1=getCov($cChr,($cSt-$cLen),($cSt-1),"q0",50,$cSample);
	$cDFL2=getCov($cChr,($cEn+1),($cEn+$cLen),"q0",50,$cSample);
	$cDFL=avgIfGr0($cDFL1,$cDFL2);
	$cDFL=$meanCov{$cChr} if $cDFL<0;
	
	# for internal try high MQ else all
	$cDIN=getCov($cChr,$cSt,$cEn,"q0",50,$cSample);
	$cDIN=getCov($cChr,$cSt,$cEn,"q0",0,$cSample) if $cDIN<0;

	if($cDIN<0){return}	
	$cDRF=($cDFL<=0)? (-1):sprintf("%.2f",$cDIN/$cDFL);
	$cDRA=($meanCov{$cChr}<=0)? (-1):sprintf("%.2f",$cDIN/$meanCov{$cChr});
	
	
	$cDRF=(-1) if $cDRF<0;
	$cDRA=(-1) if $cDRA<0;

	return;
}

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
    
#     undef @aCov; undef @aMQ;
# 	@aCov = $bwObj{$covQ}->features(-seq_id=>$qChr,-start=>$qSt,-end=>$qEn,-type=>"bin:$qBins");
# 	@aMQ = $bwObj{mq}->features(-seq_id=>$qChr,-start=>$qSt,-end=>$qEn,-type=>"bin:$qBins") if $qMQ>0;	
# 	
# # 	print STDERR "$qChr,$qSt,$qEn,$qMQ,$cSample,$qLen, qMQ: $qMQ, bins: $qBins $cSample\n";
#     
# 	for ($i=0; $i<=$#aCov;$i++){		
# 	  my ($cSt,$cEn,$cCov,$cMQ)  = ($aCov[$i]->start,$aCov[$i]->end,$aCov[$i]->mean,0);
# # 	  print STDERR "c: $cSt,$cEn,$cLen: $cCov ".($cMQ<$qMQ or $cCov<0)."\n";
# 	  
# 	  $cMQ=($aMQ[$i]->mean) if $qMQ>0;
# 	  my $cLen = ($cEn-$cSt+1); 
# # 	  print "c: $cSt,$cEn,$cLen: $cCov,$cMQ ".($cMQ<$qMQ or $cCov<0)."\n";
# 	  next if $qMQ>0 and $cMQ<$qMQ;
# 	  next if $cCov<0;	  
# 	  $sumCov+=$cCov; $fBins++;
# 	}
# # 	print "f: $fBins ".sprintf("%.1f",$sumCov/($fBins+0.0000001))."\n";
# 	
	
	#print STDERR "getCov $qChr, $qSt, $qEn, $qBins\n";
	
	my $aCov_stat = $bwObj{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	my $aMQ_stat = $bwObj{mq}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean') if $qMQ>0;	
	
	for ($i=0; $i<=$#$aCov_stat;$i++){		
		
      $cMQ=$qMQ>0 ? ${$aMQ_stat}[$i]:0;
      $cCov=${$aCov_stat}[$i];

	  next if $qMQ>0 and $cMQ<$qMQ;
	  next if $cCov<0;	  
	  $sumCov+=$cCov; $fBins++;
	  		
	}	
	

	if($fBins==0 or $fBins/$qBins<0.33333){ return (-1) }
	else{ return sprintf("%.1f",$sumCov/$fBins) }
}

sub avgIfGr0{
	my ($v1,$v2)=@_;
	$retVal=(-1);
	if($v1>0 and $v2>0){ $retVal=($v1+$v2)/2;  }
	elsif($v1>0){ $retVal=$v1;  }
	elsif($v2>0){ $retVal=$v2;  }
	return $retVal;
}
		
sub concordantSRPEs{

	my ($qChr,$qSt,$qEn,$cType,$tabix)=@_;
	$qLen=($qEn-$qSt+1);
	$qLenHalf=round($qLen/2,0);
	$FCount=0;
	my $verbose=0;
	
	$innerLimit=int($qLen*0.1);
	$innerLimit=$qLenHalf if $innerLimit<$qLenHalf;
	$innerLimit=5000 if $innerLimit>5000;
	
	$outerLimit=int($qLen*0.1);
	$outerLimit=1000 if $outerLimit<1000;
	$outerLimit=5000 if $outerLimit>5000;

	($minSt,$maxSt)=(int($qSt-$outerLimit),int($qSt+$innerLimit));
	($minEn,$maxEn)=(int($qEn-$innerLimit),int($qEn+$outerLimit));
			
	$minSt=1 if $minSt<0;
	
	print "# SRPE brID:$cBrID # $qChr,$qSt,$qEn $cType st($minSt,$maxSt) en($minEn,$maxEn) \n" if $verbose;
	

	$res = $tabix->query($qChr,$minSt,$maxSt);
	$rC=0;
	if(defined $res->get){
	while(my $line = $tabix->read($res)){
		
		# 1       10100   type=1m,2p;type1=1m;NBP1=10167;MQ1=2;chr2=2;BP2=33141456;type2=2p;NBP2=33141533;MQ2=0	
		($bChr1,$bPos1,$bName)=split("\t",$line);
		undef %t; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $t{$cTMP2[0]}=$cTMP2[1]; }
		($bChr2,$bPos2)=($t{chr2},$t{BP2});
		
		next if $qChr ne $bChr2;
		next if overlap($bPos1,$t{NBP1},$t{BP2},$t{NBP2})>20; # read 1 and read2 overlapping		
		
# 		print "$bChr1,$bPos1,$bChr2,$bPos2  $t{type1}:$t{type2} ($minEn<$bPos2 and $bPos2<$maxEn) oooooooooo \n" if $verbose;
		
		next if !($minEn<$bPos2 and $bPos2<$maxEn); 

		next if !exists($cat2type{ $t{type1}.":".$t{type2} }); ####### next if type does not exist
		$bType=$cat2type{$t{type1}.":".$t{type2}};		
		next if($cType ne $bType );		####### next if type is not satisfied

		last if $rC++>$S_SV_control_numSamples*100;
				
		print "$bChr1,$bPos1,$bChr2,$bPos2  $t{type1}:$t{type2}=$bType ($minEn<$bPos2 and $bPos2<$maxEn) ".abs($bPos1-$qSt).":".abs($bPos2-$qEn).":$qLen oooooooooo \n" if $verbose;
		
		$FCount++;	
			
	}}
	
	return $FCount;

}






sub overlap{
	my($R1,$R2,$Q1,$Q2)=@_; 
	($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	$ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	$returnVal=($ovlLen>0)? $ovlLen:0;
	return $returnVal;
}

sub median { ($a)=@_; @$a= sort { $a <=> $b } @$a; return $$a[((ceil(scalar(@$a)/2))-1)]; }



sub round{ my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec; }

sub roundA{ my ($cFreq)=@_;  if($cFreq!=1){ $RL=int(log10($cFreq)*(-1))+2; $cFreq=round($cFreq,$RL);  }  return $cFreq }
sub log10 { my $n = shift;  my $r=($n>0)? log($n)/log(10):0; return $r; }

