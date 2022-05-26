
# Author: Andre Minoche, a.minoche@garvan.org.au

# perl /g/data2/gd7/software/andmin/scripts/sv/conc_sv.pl -useVType -a /g/data2/gd7/research/NA12878/v4/SVs/joined/SV-CNV.vcf -b /g/data2/gd7/research/NA12878/v5/SVs/joined/SV-CNV.vcf 

# 	option --useVType
# 
# 	concA 
# 	concB
# 	exrtaA
# 	extraB
# 
# 	details file: coords frags, dist to breakpoint % of lenA, % of lenB
# 	simple conc
# 	multi part
# 
# 	Test: concA + extraA = total A
# 	Test: concB + extraB = total B
# 	Test: swap A B, results should be according
#   Test: twice same input -> 100% concordance perl /g/data2/gd7/software/andmin/scripts/sv/conc_sv.pl -useVType -a /g/data2/gd7/research/NA12878/v4/SVs/joined/SV-CNV.vcf -b /g/data2/gd7/research/NA12878/v4/SVs/joined/SV-CNV.vcf 
# 	input vcf or bed


use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
use Getopt::Long;
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use My::InterS qw(merge_segs);



$concDist=1000;
$cnvPercentOvl=80;
$useVType=0;



GetOptions ("useVType" => \$useVType,
			"a=s"  => \$inA,
			"b=s"  => \$inB,
			"distance=f" => \$concDist,
			"cnvPercentOvl=f" => \$cnvPercentOvl,
			"onlySTR"  => \$onlySTR, "notSTR"  => \$notSTR,
			"onlyPASS"  => \$onlyPASS, "notPASS"  => \$notPASS,
			"onlySeqDup"  => \$onlySeqDup, "notSeqDup"  => \$notSeqDup,
			"onlySmall"  => \$onlySmall, "notSmall"  => \$notSmall,
			"onlyLargeCNVs"  => \$onlyLargeCNVs,
			"highCNV"  => \$highCNV,
			"passCNV"  => \$passCNV,
			"highSV"  => \$highSV,
			"passSV"  => \$passSV
			
			
)  or die("Error in command line arguments\n");

$verbose=1;
print STDERR "Summed distance to break-point pair at which two SV are still concordant: $concDist\n";
print STDERR "A $inA\n";
print STDERR "B $inB\n";

print "A $inA\n";
print "B $inB\n";


#############################
############################# 	1. tabix file A
#############################

%allowedChr;
%chr2Max;

sub checkFilter{

	my ($X,$AorB)=@_;
	$cSample=${$$X{samples}}[0];
	
	if($onlySTR){ return 0 if !exists($$X{INFO}{CR}) or $$X{INFO}{CR}>0.2 }
	if($onlyPASS){ return 0 if !exists($$X{GTS}{FT}{$cSample}) or $$X{GTS}{FT}{$cSample} ne "PASS" }
	if($onlySeqDup){ return 0 if !exists($$X{INFO}{SEGD}) }
	if($onlySmall){ return 0 if !exists($$X{INFO}{SVLEN}) or $$X{INFO}{SVLEN}>50 }
	
	if($notSTR){ return 0 if exists($$X{INFO}{CR}) and $$X{INFO}{CR}<=0.2 }
	if($notPASS){ return 0 if exists($$X{GTS}{FT}{$cSample}) and $$X{GTS}{FT}{$cSample} eq "PASS" }
	if($notSeqDup){ return 0 if exists($$X{INFO}{SEGD}) }
	if($notSmall){ return 0 if exists($$X{INFO}{SVLEN}) and $$X{INFO}{SVLEN}<=50 }
	
	if($highCNV){ return 0 if !exists($$X{GTS}{FT}{$cSample}) or $$X{GTS}{FT}{$cSample} ne "HIGH" or $$X{INFO}{CNV}==0}# or $$X{INFO}{SVTYPE} !~ /DEL|DUP/ }
	if($passCNV){ return 0 if !exists($$X{GTS}{FT}{$cSample}) or ($$X{GTS}{FT}{$cSample} ne "PASS" and $$X{GTS}{FT}{$cSample} ne "HIGH")  or $$X{INFO}{CNV}==0}# or $$X{INFO}{SVTYPE} !~ /DEL|DUP/ }
	
	if($highSV){ return 0 if !exists($$X{GTS}{FT}{$cSample}) or $$X{GTS}{FT}{$cSample} ne "HIGH"  }
	if($passSV){ return 0 if !exists($$X{GTS}{FT}{$cSample}) or ($$X{GTS}{FT}{$cSample} ne "PASS" and $$X{GTS}{FT}{$cSample} ne "HIGH") }
	
	
	
	if($onlyLargeCNVs){ return 0 if ($$X{INFO}{SVTYPE} ne "DEL" and $$X{INFO}{SVTYPE} ne "DUP") or $$X{GTS}{FT}{$cSample} eq "LOW" or $$X{INFO}{CNV}==0 or $$X{GTS}{IDD}{$cSample}==1 or $$X{GTS}{DRA}{$cSample}==(-1) or ( $$X{INFO}{SVLEN}<100000 and exists($$X{INFO}{SEGD}) ) or ( $$X{INFO}{SVLEN}<10000 and !exists($$X{INFO}{SEGD}) )    }
# 	print ${$X->print_line};
	return 1;
}


if ($inA =~ /vcf$|vcf.gz$/){

	$I=My::VCF2->new($inA);
	
	open(INPUTA," |  sed \'s/END=/XEND=/\'  | sed \'s/=END/=XEND/\' > inputA.vcf") || die "can not open write inputA";
	print INPUTA ${$I->print_header};

	while (my $I=$I->next_line){ 
			print INPUTA ${$I->print_line};# if checkFilter($I,"A")==1;	
			next if $onlyLargeCNVs and $$V{INFO}{CNV}==0;
			($cChr,$cEnd)=($$I{CHROM},$$I{INFO}{"END"}); if($cEnd =~ /^[0-9]+$/){  $chr2Max{$cChr}+=0; $chr2Max{$cChr}=$cEnd if( $chr2Max{$cChr} < $cEnd); } # tabix gives not result if end search coord too big
	}close(INPUTA);
	print STDERR `sort_bgzip inputA.vcf; tabix -f -p vcf inputA.vcf.gz`;
	
	$AisVCF=1;$AFileEnding="vcf";
	open(IN,"gzip -d -c inputA.vcf.gz | ") || die "can not find $inA";
	while(<IN>){ next if /^#/; @_=split("\t",$_); $allowedChr{$_[0]}++; }close(IN);
	$W=My::VCF2->new("inputA.vcf.gz");
	
}elsif ($inA =~ /bed$|bed.gz$/){

	print STDERR `sort -k1,1 -k2,2n $inA | bgzip -f -c > inputA.bed.gz` if ($inA =~ /bed$/);
	print STDERR `gzip -d -c $inA  | sort -k1,1 -k2,2n | bgzip -f -c > inputA.bed.gz` if ($inA =~ /bed.gz$/);
	print STDERR `tabix -f -p bed inputA.bed.gz`;
	$AisVCF=0;$AFileEnding="bed";
	
	open(IN,"gzip -d -c $inA | ") || die "can not find $inA" if ($inA =~ /bed.gz$/);
	open(IN,"cat $inA | ") || die "can not find $inA" if ($inA =~ /bed$/);
	while(<IN>){ @_=split("\t",$_); $allowedChr{$_[0]}++; }close(IN);	
	
}

#print STDERR scalar(keys %allowedChr)." chr in input file\n";


if ($inB =~ /vcf$|vcf.gz$/){	
	$BisVCF=1; $BFileEnding="vcf";
	$V=My::VCF2->new($inB);
}elsif ($inB =~ /bed$|bed.gz$/){
	$BisVCF=0; $BFileEnding="bed";	
}

$tabix = new Tabix(-data => "inputA.$AFileEnding.gz",-index => "inputA.$AFileEnding.gz.tbi");



open(DET,">details.txt") || die "can not open details.txt";
open(CONCA,">CONC_A.$AFileEnding") || die "can not open details.$AFileEnding";
open(CONCB,">CONC_B.$BFileEnding") || die "can not open details.$BFileEnding";
open(EXTRAA,">EXTRA_A.$AFileEnding") || die "can not open details.$AFileEnding";
open(EXTRAB,">EXTRA_B.$BFileEnding") || die "can not open details.$BFileEnding";

print CONCB ${$V->print_header} if $BisVCF;
print EXTRAB ${$V->print_header} if $BisVCF;
print CONCA ${$W->print_header} if $AisVCF;
print EXTRAA ${$W->print_header} if $AisVCF;


#############################
############################# 2. loop through B
#############################

our %Aconc=();
if($BisVCF){
    
	while (my $V=$V->next_line){
		
		next if checkFilter($V,"B")==0;
		
		$allVariantsinA++;
	    my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType)=obtainCoordsFromVCF($V);
	    push @{$stat{countB}}, $cLen;
# 		print STDERR ${$V->print_line} if $verbose;
		
		#if(!exists($allowedChr{$cChr1}) ){print EXTRAB ${$V->print_line}; next}
		
		$passCat="";
		$resConc=conCheck($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType);
		
		###### print conc or extra
		if ($resConc eq "conc"){
			push @{$stat{concB}}, $cLen;
			print CONCB ${$V->print_line};	
		}else{
			push @{$stat{extraB}}, $cLen;	
			
			if ($cChr1 ne $cChr2){
	    		print DET "extra\tB\t$passCat\t$cChr1:$cPos1,$cChr2:$cPos2,$aType\t\t".${$V->print_line}."";
	    	}else{
	    		print DET "extra\tB\t$passCat\t$cChr1:$cPos1-$cPos2,$aType\t\t".${$V->print_line}."";
	    	}print EXTRAB ${$V->print_line};
		}	
	}
}else{ # for bed

	open(IN,"gzip -d -c $inB | ") || die "can not find $inA" if ($inB =~ /bed.gz$/);
	open(IN,"cat $inB | ") || die "can not find $inA" if ($inB =~ /bed$/);
	while (<IN>){
	
		chomp; my ($cChr1,$cPos1,$cPos2,$cType)=split("\t",$_);
# 		print "$_\n";
		$cPos1++;
		$cChr2=$cChr1;
		$cLen=$cPos2-$cPos1+1;
		push @{$stat{countB}}, $cLen;
		
		if(!exists($allowedChr{$cChr1}) ){print EXTRAB $_."\n"; next}
		$passCat="";
		$resConc=conCheck($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType);
		
		###### print conc or extra
		if ($resConc eq "conc"){
			print CONCB $_."\n";
			push @{$stat{concB}}, $cLen;	
		}else{
			print EXTRAB $_."\n";
						
			if ($cChr1 ne $cChr2){
	    		print DET "extra\tB\t$passCat\t$cChr1:$cPos1,$cChr2:$cPos2,$aType\t\t\n";
	    	}else{
	    		print DET "extra\tB\t$passCat\t$cChr1:$cPos1-$cPos2,$aType\t\t\n";
	    	}
			push @{$stat{extraB}}, $cLen;	
		}
	}
	close(IN);

}

#############################
############################# 3. loop through A
#############################
if($AisVCF){
	while (my $W=$W->next_line){
		
		$allVariantsinA++;
	    my ($aChr1,$aPos1,$aChr2,$aPos2,$aLen,$aType)=obtainCoordsFromVCF($W);
		push @{$stat{countA}}, $aLen;
		
	    if (  exists($Aconc{"$aChr1,$aPos1,$aChr2,$aPos2,$aType"})  ){
	    	print CONCA ${$W->print_line};
	    	push @{$stat{concA}}, $aLen;	
	    }else{
	    	print EXTRAA ${$W->print_line};
	    	push @{$stat{extraA}}, $aLen;
	    	if ($aChr1 ne $aChr2){
	    		print DET "extra\tA\t\t$aChr1:$aPos1,$aChr2:$aPos2,$aType\t\t\n";
	    	}else{
	    		print DET "extra\tA\t\t$aChr1:$aPos1-$aPos2,$aType\t\t\n";
	    	}
	    }
	}
}else{
	open(IN,"gzip -d $inA | ") || die "can not find $inA" if ($inA =~ /bed.gz$/);
	open(IN,"<$inA") || die "can not find $inA" if ($inA =~ /bed$/);
	while (<IN>){
	
		chomp; my ($aChr1,$aPos1,$aPos2,$aType)=split("\t",$_);
		$aPos1++; $aLen=($aPos2-$aPos1+1); $aChr2=$aChr1;
		push @{$stat{countA}}, $aLen;
		
		###### print conc or extra
		if (  exists($Aconc{"$aChr1,$aPos1,$aChr2,$aPos2,$aType"})  ){
	    	print CONCA $_."\n";
	    	push @{$stat{concA}}, $aLen;	
	    }else{
	    	print EXTRAA $_."\n";
	    	push @{$stat{extraA}}, $aLen;
	    	if ($aChr1 ne $aChr2){
	    		print DET "extra\tA\t\t$aChr1:$aPos1,$aChr2:$aPos2,$aType\t\t\n" ;
	    	}else{
	    		print DET "extra\tA\t\t$aChr1:$aPos1-$aPos2,$aType\t\t\n" ;
	    	}
	    }
	    
	}
	close(IN);
}



## calculate current coverage

#############################
############################# 4. print conc summary
#############################


# do a histogram with givin intervals
@categories=("countA","countB","concA","concB","extraA","extraB");
@intervals=(0,500,1000,5000,10000,50000,100000,500000,5000000,100000000,250000000);

print join("\t",("fragmentSize",@categories))."\n";

for ($i=1;$i<=$#intervals;$i++){

	($l1,$l2)=($intervals[$i-1],$intervals[$i]);

	undef @tmp2;
	foreach $cCat (@categories){
	 @tmp=grep { $_ > $l1 && $_ <= $l2} @{$stat{$cCat}};
	 push @tmp2, scalar(@tmp);
	}
	print  join("\t",(formatNumber($l1)."-".formatNumber($l2),@tmp2))."\n";
	
}
undef @tmp2; foreach $cCat (@categories){ push @tmp2, scalar(@{$stat{$cCat}}); }
print  join("\t",("Total",@tmp2))."\n";

#undef @tmp2; foreach $cCat (@categories){ print join("\t",("$cCat: ",grep { $_ < $intervals[0] || $_ > $intervals[$#intervals]} @{$stat{$cCat}}))."\n"; }


sub formatNumber {
	($n)=@_;
	$n=~ s/000000000$/G/;
	$n=~ s/000000$/M/;
	$n=~ s/000$/k/;
	return $n;
}

# fragmentSize	countA	countB	AconcToB	BconcToA	Aextra	Bextra	Sensitivity	Precision
# 1-500	246	2848	41	19	205	2829	1%	17%
# 500-1k	127	453	34	28	93	425	6%	27%
# 1k-5k	238	526	71	60	167	466	11%	30%
# 5k-10k	64	154	19	19	45	135	12%	30%
# 10k-50k	103	38	15	16	88	22	42%	15%
# 50k-100k	28	4	3	1	25	3	25%	11%
# 100k-500k	60	8	2	4	58	4	50%	3%
# 500k-5M	21	0	0	0	21	0	#DIV/0!	0%
# 5M-200M	10	0	0	0	10	0	#DIV/0!	0%
# >200M	0	0	0	0	0	0	#DIV/0!	#DIV/0!
# 	897	4031	185	147	712	3884		





#############################
############################# Functions
#############################

sub	conCheck {
	
	my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType)=@_;
	
	
# 	print "--- $cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType\n" if $verbose;

	####### 1 compare distance to breakpoint
	@testCoords=([$cChr1,$cPos1,$cChr2,$cPos2,$cType]);#,[$cChr2,$cPos2,$cChr1,$cPos1,$cType]);
	$concFound=0;
	
	foreach (@testCoords){
	
		($BChr1,$BPos1,$BChr2,$BPos2,$BType)=@$_;
		$cPos1P=($BPos1+$concDist);
		$cPos1M=($BPos1-$concDist);		
# 		print DET "    ---test $BChr1,$cPos1M,$cPos1P\n" if $verbose;
		$res = $tabix->query($BChr1,$cPos1M,$cPos1P);
		
		$passCat="tabix";
		if(defined $res->get){
		while(my $line = $tabix->read($res)){ #1       869477  870222  name=Lumpy_28;sample=FR05812673;SVTYPE=DEL;SVLEN=745;TOOL=Lumpy;SR=0;PE=8;DRF=0.12;IDD=.;       0       +       869477  870222  DEL     1       746     0
		
# 			print DET "      tabix1 --- $line\n" if $verbose;
			
			if(!$AisVCF){chomp($line); ($AChr1,$APos1,$APos2,$AType)=split("\t",$line); $APos1++; $ALen=($APos2-$APos1+1); $AChr2=$AChr1; }
			else{  $W->parse_line(\$line); ($AChr1,$APos1,$AChr2,$APos2,$ALen,$AType,$formatAll,$gtAll)=obtainCoordsFromVCF($W);    }
			
			$passCat="type"; next if $BType ne $AType and $useVType; # type DEL, INV, .. has to match
			$passCat="chr"; next if $AChr1 ne $BChr1 or $AChr2 ne $BChr2; # chr have to match
		
			$cDist1=abs($APos1-$BPos1); $cDist2=abs($APos2-$BPos2); $cSumDist=($cDist1+$cDist2);
			
			$passCat="distc.$cDist1.$cDist2"; next if $cDist1>$concDist or $cDist2>$concDist; # pos not further away than concDist
			$passCat="ovl"; next if $AChr1 eq $AChr2 and overlap($APos1,$APos2,$BPos1,$BPos2)==0; # of on same chromosome, there must be a overlap. This is important for the small events < concDist
		
			#next if exists($Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"});
			$Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"}++;
			$Bconc{"$BChr1,$BPos1,$BChr2,$BPos2,$BType"}++;
			$passCat="conc";
			$concFound++;
			($cLen,$ALen)=(0,0) if $BChr1 ne $BChr2;
			if($cType eq "BND"){
				print DET "concordant\tB\tdistance\t$BChr1:$BPos1,$BChr2:$BPos2,$BType\tD:$cSumDist\t$formatAll\t$gtAll\n" if $Bconc{"$BChr1,$BPos1,$BChr2,$BPos2,$BType"}==1;
				print DET "concordant\tA\tdistance\t$AChr1:$APos1,$AChr2:$APos2,$AType\tD:$cSumDist\t$formatAll\t$gtAll\n" if $Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"}==1;
				print DET "det_concA\tdistance\t$cSumDist\tA:$AChr1:$APos1,$AChr2:$APos2,$AType\tB:$BChr1:$BPos1,$BChr2:$BPos2,$BType\t$cDist1,$cDist2\t$cLen,$ALen\n";
			}else{
				print DET "concordant\tB\tdistance\t$BChr1:$BPos1-$BPos2,$BType\tD:$cSumDist\t$formatAll\t$gtAll\n" if $Bconc{"$BChr1,$BPos1,$BChr2,$BPos2,$BType"}==1;
				print DET "concordant\tA\tdistance\t$AChr1:$APos1-$APos2,$AType\tD:$cSumDist\t$formatAll\t$gtAll\n" if $Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"}==1;
				print DET "det_concA\tdistance\t$cSumDist\tA:$AChr1:$APos1-$APos2,$AType\tB:$BChr1:$BPos1-$BPos2,$BType\t$cDist1,$cDist2\t$cLen,$ALen\n";
			}
		}
		}
# 		last if $concFound;
	}
# 	print $concFound." dist \n";
	return "conc" if $concFound>0;
		
	####### 2 for CNVs
	
	if(($cType =~ /DUP|DEL/ and $useVType) or !$useVType){
		$concFound=0;
		
		
		($BChr1,$BPos1,$BChr2,$BPos2,$BLen,$BType)=($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType);
		#print "    ---x2 $BChr1,$BPos1,$BPos2 $cType\n" if $verbose;
		
		$min=($BPos1-$BLen*100);$min=0 if $min<=0;
		$max=($BPos2+$BLen*100); $max=$chr2Max{$BChr1} if exists($chr2Max{$BChr1}) and $max>$chr2Max{$BChr1};
		$res = $tabix->query($BChr1,$min,$max);
		
		$passCat="tabix2";
		undef @allOvl; # start2end
		
		if(defined $res->get){
		while(my $line = $tabix->read($res)){ #1       869477  870222  name=Lumpy_28;sample=FR05812673;SVTYPE=DEL;SVLEN=745;TOOL=Lumpy;SR=0;PE=8;DRF=0.12;IDD=.;       0       +       869477  870222  DEL     1       746     0
		
			if(!$AisVCF){chomp($line); ($AChr1,$APos1,$APos2,$AType)=split("\t",$line); $APos1++; $ALen=($APos2-$APos1+1); $AChr2=$AChr1; }
			else{  $W->parse_line(\$line); ($AChr1,$APos1,$AChr2,$APos2,$ALen,$AType,$formatAll,$gtAll)=obtainCoordsFromVCF($W);   }
# 			print "      tabix2 --- $AChr1,$APos1,$AChr2,$APos2,$ALen,$AType\n";
			
			$passCat="tabix_type";
			next if $BType ne $AType and $useVType; # type DEL, INV, .. has to match
			
			$passCat="tabix_chr";
			next if $AChr1 ne $BChr1 or $AChr2 ne $BChr2; # chr have to match
			
			$passCat="tabix_ovl";
			$cOvl=overlap($APos1,$APos2,$BPos1,$BPos2);
			next if $cOvl==0;
			
			$BOvlPc=round($cOvl/$BLen*100,0);
			$AOvlPc=round($cOvl/$ALen*100,0);
			$cSumPc=($BOvlPc+$AOvlPc);
# 			print DET "$cOvl/$BLen   $cOvl/$ALen $APos1,$APos2,$BPos1,$BPos2 \n";
 
			$passCat="ovl2";

			push @allOvl, [$AChr1,$APos1,$APos2,"$AChr1,$APos1,$AChr2,$APos2,$AType,$BOvlPc,$AOvlPc,$cSumPc"];

# 			print "$BOvlPc $AOvlPc $cnvPercentOvl -- pass\n";				
			if($AOvlPc>=$cnvPercentOvl & $BOvlPc>=$cnvPercentOvl){
				$Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"}++;
				
				print DET "concordant\tA\toverlap\t$AChr1:$APos1-$APos2,$AType\tA:$AOvlPc,B:$BOvlPc\t$formatAll\t$gtAll\n" if $Aconc{"$AChr1,$APos1,$AChr2,$APos2,$AType"}==1;
				print DET "det_concA\toverlap\tsumA:$AOvlPc\tA:$AChr1:$APos1-$APos2,$AType\tB:$BChr1:$BPos1-$BPos2,$BType\tA:$AOvlPc,B:$BOvlPc\tA:$ALen,B:$cLen\n";
				
				
			}
		}	
		}
		
		if(@allOvl>0){
		
			merge_segs(\@allOvl); ## merge if possible, repeat until last merge did not decreased array size
			$BPcOvl=calcOvlLen(\@allOvl,$BPos1,$BPos2);
			
			if(scalar(@allOvl)==1){	
				$AOvlPc=round(overlap($allOvl[0][1],$allOvl[0][2],$BPos1,$BPos2)/($allOvl[0][2]-$allOvl[0][1]+1)*100,0);
# 				print DET $allOvl[0][3]."\n";
			}else{ $AOvlPc="M"; }

			
			if($BPcOvl >= $cnvPercentOvl){
			
				$cDiscOrConc="det_concB";
				$concFound++;	
				$Bconc{"$BChr1,$BPos1,$BChr2,$BPos2,$BType"}++;
				print DET "concordant\tB\toverlap\t$BChr1:$BPos1-$BPos2,$BType\tA:$AOvlPc,B:$BPcOvl\n" if $Bconc{"$BChr1,$BPos1,$BChr2,$BPos2,$BType"}==1;
			
			}else{	
				$cDiscOrConc="det_discB";
			}
						
			# assign concordant to A
			$oC=0;
			foreach $cArr (@allOvl){ # foreach segment
				foreach $cStr (@{$$cArr[3]}){ # foreach merged segment
					($AChr1,$APos1,$AChr2,$APos2,$AType,$BOvlPc,$AOvlPc,$cSumPc)=split(",",$cStr);
					$oC++;
					print DET "$cDiscOrConc\toverlapSum\tm$oC,sumB:$BPcOvl\tA:$AChr1:$APos1-$APos2,$AType\tB:$BChr1:$BPos1-$BPos2,$BType\tA:$AOvlPc,B:$BOvlPc\tA:$ALen,B:$cLen\tA:$formatAll\tA:$gtAll\n";
				}
			}
			

		}
		
		
# 		print $concFound." ovl \n";
		return "conc" if $concFound>0;	
	}
	
	return "extra";
}


sub obtainCoordsFromVCF{
	
	($X)=@_;
	($aChr1,$aPos1)=($$X{CHROM},$$X{POS});
	
	if ($$X{"ALT"} =~ /del|DEL/ or $$X{"INFO"}{"SVTYPE"} =~ /del|DEL/ ){$aType="DEL"}
	elsif ($$X{"ALT"} =~ /dup|DUP/ or $$X{"INFO"}{"SVTYPE"} =~ /dup|DUP/ ){$aType="DUP"}
	elsif ($$X{"ALT"} =~ /inv|INV/ or $$X{"INFO"}{"SVTYPE"} =~ /inv|INV/ ){$aType="INV"}
	else{$aType="BND"}
	
	if ($$X{"INFO"}{"SVTYPE"} eq "BND"){
		# 		next if exists($$X{"INFO"}{"SECONDARY"});
		if($$X{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){

			($aChr2,$aPos2)=($1,$2);
			if($aChr1 eq $aChr2){ 	
				$$X{"INFO"}{"SVLEN"}=abs($aPos1-$aPos2)+1;
				$$X{"INFO"}{"XEND"}=$aPos2;
			}
		}else{die "can not interprete ALT field".$$X{"all"} }
	}else{
		$aPos2=(exists($$X{"INFO"}{"XEND"}))? $$X{"INFO"}{"XEND"}:$$X{"INFO"}{"END"};
		($aChr2)=($$X{CHROM});
	}
	$aLen=abs($aPos2-$aPos1)+1;

	$$X{"INFO"}{"SVLEN"}=$aLen if !exists($$X{"INFO"}{"SVLEN"}) and $aChr1 eq $aChr2;
	$$X{"INFO"}{"SVLEN"}=abs($$X{"INFO"}{"SVLEN"}) if exists($$X{"INFO"}{"SVLEN"});
# 	print  "$aChr1,$aPos1,$aChr2,$aPos2,$aLen,$aType oo\n";
	return((($aChr1,$aPos1,$aChr2,$aPos2,$aLen,$aType,$$X{"FORMAT_all"},$$X{"GT_all"})));
}



sub overlap{ my($R1,$R2,$Q1,$Q2)=@_; ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;}





sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}

sub calcOvlLen{
	
	($arr,$rSt,$rEn)=@_;

	my $sumOvl=0;
	
	for ($i=0; $i<=$#$arr;$i++){
		$sumOvl+=overlap($$arr[$i][1],$$arr[$i][2],$rSt,$rEn);			
	}
# 	print STDERR " $sumOvl/($rEn-$rSt+1)\n";
	return round($sumOvl/($rEn-$rSt+1)*100,0);
}





























