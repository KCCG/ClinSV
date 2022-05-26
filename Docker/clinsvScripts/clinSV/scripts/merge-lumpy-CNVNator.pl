
# Author: Andre Minoche, a.minoche@garvan.org.au


use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use List::BinarySearch::XS qw( binsearch_pos );
use Sort::Naturally;
use My::VCF2; 

$CNVnatorStem=$ARGV[0];
$samples=$ARGV[1];
$lumpyVCF=$ARGV[2];
$OutStemMerged=$ARGV[3];
@Asamples=split(",",$samples);


$ovlCutOff=0.7;
$ovlCutOffLumpy=0.7;
$ovlTolerance=200;
$ovlToleranceLumpy=100;


$verbose=0;

############################## merge CNVnator calls of multiple samples

# the frequency must be the same
# 2       89917200        89923300        0.63_6k 2.61266e-11     0.63    0       6100

%CNVmaster={};
%CNVslave={};
%iCNVmaster={};
$inFileCount=0;

foreach $cSample (@Asamples){
	
	$inFileCount++;
	$cinputRows=0;
	$cInF=$CNVnatorStem; $cInF=~ s/XX_SAMPLEID_XX/$cSample/g;
	
	print STDERR "reading $cSample \n";
	open(IN1, "<$cInF") || die "$cInF nicht gefunden";
	while(<IN1>){ chomp;
		
		next if length($_)==0;
		$cinputRows++;
		($cChr,$cSt,$cEn,$cName)=split("\t",$_);
		undef %h; foreach $cTMP (split(";",$cName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=$cTMP2[1]; }
		($cPval,$cFreq,$cLen,$cFT)=($h{"PVAL"},$h{"FREQ"},$h{"SVLEN"},$h{"FT"});
		
		$cSt++; if($cSt>$cEn){  ($cSt,$cEn)=($cEn,$cSt); }
		
		
		$notMerged=tryToMerge($cChr,$cSample,$cSt,$cEn,$cLen,$cFreq,$cPval,$cName,$cFT);
		
		push @{$CNVslave{$cChr}}, [$cSt,$cEn,$cLen,$cFreq,[$cSample],[$cSt],[$cEn],[$cLen],[$cFreq],[$cPval],[$cName],[$cFT]] if $notMerged;
		$lastElem=$#{$CNVslave{$cChr}};
		print "0:  $cFreq ${$CNVslave{$cChr}}[$lastElem][3] $cChr:$cSt-$cEn\n" if $verbose;
		
	}close(IN1);
		
	$rowsMaster=mergeMaster();
	print  STDERR "  $cSample $cinputRows rows => $rowsMaster cluster \n";

}

# dumpMaster();
sub dumpMaster{
	#print master table
	foreach $cChr (nsort keys %CNVmaster){
	
	# 	print join("\n",@{$iCNVmaster{$cChr}})."\n";
		$aChr=$CNVmaster{$cChr};
		for ($i=0; $i<=$#$aChr;$i++){
			$cL=$$aChr[$i];
# 			next if exists($CNVMNatorLineBlackList{$cChr}{$i});
			@closToPrint=($cChr,@$cL[0..3]);
			push @closToPrint,join(",",@{$$cL[$_]}) for (4..$#$cL);
			print OUT join("\t",@closToPrint)."\n";
		}
	}
}

sub mergeMaster{
	# merge to master table
	foreach $cChr (keys %CNVslave){ push @{$CNVmaster{$cChr}}, @{$CNVslave{$cChr}}; } %CNVslave={};
	
	# sort master table
	foreach $aChr (values %CNVmaster){ @$aChr=sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @$aChr; }
	
	# create index for master
	foreach $cChr (keys %CNVmaster){@{$iCNVmaster{$cChr}}=(); $iChr=$iCNVmaster{$cChr}; foreach $cL (@{$CNVmaster{$cChr}}){ push @$iChr,$$cL[0]; } }
	
	# count rows master
	$rowsMaster=0; foreach $cChr (keys %CNVmaster){ $rowsMaster+=scalar(@{$CNVmaster{$cChr}}); }
	
	return $rowsMaster;
}

sub tryToMerge{

	($cChr,$cSample,$cSt,$cEn,$cLen,$cFreq,$cPval,$cName,$cFT)=@_;
	return 1 if !exists($iCNVmaster{$cChr});
	
	($iSt,$iEn)=binsearch_range(\@{$iCNVmaster{$cChr}},($cSt-$cLen),$cEn);
	return 1 if $iSt<0; # no overlap
	
	$aChr=$CNVmaster{$cChr}; # pointer to current chromosome array 
# 	print "$cSt,$cEn : $iSt,$iEn : ".scalar(@$aChr)."\n";
	for my $i ($iSt..$iEn){
		
		$cL=$$aChr[$i]; # pointer to current line
		($cSt2,$cEn2,$cLen2,$cFreq2)=@$cL[0..3];
		$cOvlLen=overlap($cSt,$cEn,$cSt2,$cEn2);
# 		print "$cSt,$cEn,$cLen,$cFreq  $cSt2,$cEn2,$cLen2,$cFreq2 \n";
		if( (($cOvlLen/$cLen) >= $ovlCutOff or abs($cOvlLen-$cLen)<=$ovlTolerance)
		and (($cOvlLen/$cLen2) >= $ovlCutOff or abs($cOvlLen-$cLen2)<=$ovlTolerance)  
		and (($cFreq>=1 and $cFreq2>=1) or ($cFreq<1 and $cFreq2<1))    ){
# 		if( ($cOvlLen/$cLen) >= $ovlCutOff  and ($cOvlLen/$cLen2) >= $ovlCutOff ){
# 			print  "mergeToCluster $i:  ($cOvlLen/$cLen) >= $ovlCutOff  and ($cOvlLen/$cLen2) >= $ovlCutOff and (($cFreq>=1 and $cFreq2>=1) or ($cFreq<1 and $cFreq2<1)) \n";
			mergeLine($cL,$cSample,$cSt,$cEn,$cLen,$cFreq,$cPval,$cName,$cFT);
			
			return 0;
		}
	}
	return 1;
}
	
sub mergeLine{

	my ($cL,$cSample,$cSt,$cEn,$cLen,$cFreq,$cPval,$cName,$cFT)=@_;
	
	### [$cSt,$cEn,$cLen,$cFreq,[$cSample],[$cSt],[$cEn],[$cLen],[$cFreq],[$cPval]]
	
	$sF=4; foreach $cEl (($cSample,$cSt,$cEn,$cLen,$cFreq,$cPval,$cName,$cFT)){ push @{$$cL[$sF]},$cEl; $sF++; }
	
	$$cL[0]=sprintf("%.0f",average(\@{$$cL[5]})); # new average start
	$$cL[1]=sprintf("%.0f",average(\@{$$cL[6]})); # ... end
	$$cL[2]=sprintf("%.0f",($$cL[1]-$$cL[0]+1));  # ... len
	$$cL[3]=average(\@{$$cL[8]}); # ... freq
	
	
	
}



############################## go through lumpy vcf and merge with CNVnator calls 
$pos2sample={};
%CNVMNatorLineBlackList;
print STDERR "parse lumpy VCF...\n";
open(OUT, ">$OutStemMerged.vcf") || die "$OutStemMerged.vcf nicht gefunden";

if(! -e $lumpyVCF ){ die "can not find vcf file\n"; }

my $VCFObj=My::VCF2->new($lumpyVCF);

$VCFObj->add_header("##FORMAT=<ID=CNRD,Number=1,Type=Float,Description=\"CNVnator read depth normalized to 1 \">");
$VCFObj->add_header("##FORMAT=<ID=CNP,Number=1,Type=Float,Description=\"CNVnator e-val by t-test\">");

$VCFObj->add_header("##INFO=<ID=TOOL,Number=1,Type=String,Description=\"Varaint detection tool: Lumpy and/or CNVnator\">");
$VCFObj->add_header("##INFO=<ID=EVTYPE,Number=1,Type=String,Description=\"Type of Lumpy evidence contributing to the variant call\">");
# $VCFObj->add_header("##INFO=<ID=CNR,Number=1,Type=Float,Description=\"CNVnator copy number ratio normalized to 1\">");

$VCFObj->add_header("##source=ClinSV");


$VCFObj->remove_header("FORMAT","CNRD");

$VCFObj->remove_header("FORMAT","IDD");
$VCFObj->remove_header("FORMAT","ICN");
$VCFObj->remove_header("FORMAT","BD");
$VCFObj->remove_header("FORMAT","DRFN");
$VCFObj->remove_header("INFO","CPE");
$VCFObj->remove_header("INFO","CSR");
$VCFObj->remove_header("INFO","OL");



# $VCFObj->add_header("##INFO=<ID=Q0,Number=1,Type=String,Description=\" % percentage of positions with MQ=0 within CNV (form CNVnator)\">");
			
print OUT ${$VCFObj->print_header};

while (my $V=$VCFObj->next_line){ # loop through vcf file
  	
  	($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{"INFO"}{"END"});
  	
  	$aSamples=\@{$$V{samples}};
  	if($$V{"INFO"}{"SVTYPE"} eq "BND" or ($$V{"INFO"}{"CNV"} == 0 and $$V{"INFO"}{"SVTYPE"} eq "DEL") ){ # BNDs and deletion copy neutral events need at least for lines of evidence 	
  		$PASS=0; foreach $cSample (@{$$V{samples}}){ $PASS=1 if $$V{GTS}{SU}{$cSample}>=4;  }
  		next if $PASS==0;
  		foreach $cSample (@{$$V{samples}}){ $$V{GTS}{GT}{$cSample}=($$V{GTS}{SU}{$cSample}>=4)? "1/1":"0/0";  }
  	 }else{
  	 	$PASS=0; foreach $cSample (@{$$V{samples}}){ $PASS=1 if exists($$V{GTS}{FT}) and $$V{GTS}{FT}{$cSample} eq "PASS";  }
  		next if $PASS==0;
  		foreach $cSample (@{$$V{samples}}){ $$V{GTS}{GT}{$cSample}=($$V{GTS}{FT}{$cSample} eq "PASS")? "1/1":"0/0";  }
  	 }
  	
  	#adjust genotype to het 0/1 if CNEUTR
  	if ($$V{"INFO"}{"CNV"} == 1){
  		foreach $cSample (@{$$V{samples}}){ 
			if($$V{"INFO"}{"SVTYPE"} eq "DEL" and $$V{GTS}{DRA}{$cSample} >=0.2 and $$V{GTS}{DRA}{$cSample} <0.8){
				$$V{GTS}{GT}{$cSample}="0/1"; 
			}elsif($$V{"INFO"}{"SVTYPE"} eq "DUP" and $$V{GTS}{DRA}{$cSample} >1.2 and $$V{GTS}{DRA}{$cSample} <=1.75){
				$$V{GTS}{GT}{$cSample}="0/1"; 
			}
  		}
	}
	
  	
  	
	$$V{"INFO"}{"TOOL"}="Lumpy";
	$$V{"ID"}="Lumpy_".$$V{"ID"};
	$$V{"INFO"}{"EVTYPE"}.="PE".$$V{"INFO"}{"PE"}."," if $$V{"INFO"}{"PE"}>0;
	$$V{"INFO"}{"EVTYPE"}.="SR".$$V{"INFO"}{"SR"}."," if $$V{"INFO"}{"SR"}>0;
	$$V{"INFO"}{"EVTYPE"}.="RD,";
	
	if(($$V{"INFO"}{"SVTYPE"} eq "DEL" or $$V{"INFO"}{"SVTYPE"} eq "DUP" ) and $$V{"INFO"}{"CNV"} == 1){ # do not merge copy neutral events del or dupes, since they are not supposed to be called by CNVnator

		($cChr,$cSt,$cEn,$SVTYPE,$cLen)=($$V{"CHROM"},$$V{"POS"},$$V{"INFO"}{"END"},$$V{"INFO"}{"SVTYPE"},abs($$V{"INFO"}{"SVLEN"}));
		if($cSt>$cEn){  ($cSt,$cEn)=($cEn,$cSt); }
# 		print "$cChr,$cSt,$cEn,$SVTYPE,$cLen\n";
		$matchLine=ovlCNVnatorLine($cChr,$cSt,$cEn,$SVTYPE,$cLen);
	
		if($matchLine>=0){
			$CNVMNatorLineBlackList{$cChr}{$matchLine}++;
			$lumpyLinesMerged++;
			$$V{"INFO"}{"EVTYPE"}=~ s/RD,/RD2,/g;
			$$V{"INFO"}{"TOOL"}.=",CNVnator";
			
			#add the CNVnator fields to current vcf line
			$aChr=$CNVmaster{$cChr};
			$cL=$$aChr[$matchLine];
			undef %tmpF;
			
			# for each samples
			for ($k=0; $k<=$#{$$cL[4]};$k++){ # [$cSt,$cEn,$cLen,$cFreq,[$cSample],[$cSt],[$cEn],[$cLen],[$cFreq],[$cPval],[$cName],[$cFT]]
				($cSample2,$cFreq2,$cPval2,$cName,$cFT)=(${$$cL[4]}[$k],${$$cL[8]}[$k],${$$cL[9]}[$k],${$$cL[10]}[$k],${$$cL[11]}[$k]);					
				$$V{GTS}{CNRD}{$cSample2}=$cFreq2;
				$$V{GTS}{CNP}{$cSample2}=$cPval2;
				
				# if lumpy not PASS, let PASS by CNVnator if CNVnator is PASS
				if ($$V{GTS}{FT}{$cSample2} ne "PASS" and $cFT eq "PASS"){
				$$V{GTS}{FT}{$cSample2}="PASS";
				$$V{GTS}{GT}{$cSample2}="1/1";
					#adjust genotype to het 0/1 if CNV
					if($$V{"INFO"}{"SVTYPE"} eq "DEL" and $$V{GTS}{DRA}{$cSample2} >=0.2 and $$V{GTS}{DRA}{$cSample2} <0.8){
						$$V{GTS}{GT}{$cSample2}="0/1"; 
					}elsif($$V{"INFO"}{"SVTYPE"} eq "DUP" and $$V{GTS}{DRA}{$cSample2} >1.2 and $$V{GTS}{DRA}{$cSample2} <=1.75){
						$$V{GTS}{GT}{$cSample2}="0/1"; 
					}
				}
			}


		}
	}	
	
	chop($$V{"INFO"}{"EVTYPE"});
	print OUT ${$V->print_line};
	
	print STDERR "$c lines processed           \r" if ($c++)%100000==0;
}


printUnmergedCNVN();
close(OUT);
print STDERR "  $lumpyLinesMerged lumpy calls merged with merged CNVnator calls\n";

sub printUnmergedCNVN{
	#print master table
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  FR05812618      FR05812626      FR05812641      FR05812650      FR05812658      FR05812673
	#1       963914  60      N       <DEL>   .       .       SVTYPE=DEL;SVLEN=-545;END=964459;STRANDS=+-:56;IMPRECISE;CIPOS=-7,61;CIEND=-77,10;CIPOS95=-4,7;CIEND95=-14,2;SU=56;PE=56;SR=0;LOCATION=1:963914-964459;LEN=545  GT:SU:PE:SR:DRF:DRS:DBP:DIN:IDD:PASS:CNF:CNP    ./.:0:0:0:0.48:1.41:29:13:.:0:.:.       ./.:9:9:0:0.20:0.42:24:5:.:1:.:.        ./.:21:21:0:0.00:0.00:26:0:.:1:0.00:0.129702    ./.:0:0:0:0.34:0.71:30:10:.:0:.:.       ./.:10:10:0:0.18:0.37:28:5:.:1:.:.      ./.:16:16:0:0.00:0.00:20:0:.:1:0.00:0.152488
	@outFTypes=("GT","CNP","SU","PAFSU","PE","SR","DBP","DRF","DRA","MQBP","FT"); # exclude "CNRD"

	foreach $cChr (nsort keys %CNVmaster){
	
	# 	print join("\n",@{$iCNVmaster{$cChr}})."\n";
		$aChr=$CNVmaster{$cChr};
		for ($i=0; $i<=$#$aChr;$i++){
			$cL=$$aChr[$i];
			next if exists($CNVMNatorLineBlackList{$cChr}{$i});
			
# 			@closToPrint=($cChr,@$cL[0..3]);
# 			push @closToPrint,join(",",@{$$cL[$_]}) for (4..$#$cL);
# 			print OUT join("\t",@closToPrint)."\n";
			undef %tmpF; $avgPVal=0;$avgPValC=0; undef @fredS; undef @MQS; undef @PCSD; undef @CPE;undef @CSR;

			for ($k=0; $k<=$#{$$cL[4]};$k++){ # [$cSt,$cEn,$cLen,$cFreq,[$cSample],[$cSt],[$cEn],[$cLen],[$cFreq],[$cPval]]
				($cSample,$cFreq,$cPval,$cName)=(${$$cL[4]}[$k],${$$cL[8]}[$k],${$$cL[9]}[$k],${$$cL[10]}[$k]);	
				$tmpF{$cSample}{"freq"}=$cFreq;
				$tmpF{$cSample}{"cPval"}=$cPval;
				$avgPValC++;$avgPVal+=$cPval;
				push @fredS, fredQ($cPval);
				
				#Q0=0.000463271;DBP=15;DRF=0.4;FT=PASS;GC=78;MQ=37.2;PE=0;MQBP=60;SR=0;
				undef %h; foreach $cTMP (split(";",$cName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=$cTMP2[1]; }
				push @MQS, $h{MQ};
				push @PCSD, $h{PCSD};
				push @ACPE, $h{CPE};
				push @ACSR, $h{CSR};
				$tmpF{$cSample}{"PE"}=$h{PE}; # @outFTypes=("GT","CNRD","CNP","SU","PAFSU","PE","SR","DBP","DRF","DRA","MQBP","FT");
				$tmpF{$cSample}{"SR"}=$h{SR};
				$tmpF{$cSample}{"SU"}=$h{SU};
				$tmpF{$cSample}{"PAFSU"}=$h{PAFSU};
				$tmpF{$cSample}{"DBP"}=$h{DBP};
				$tmpF{$cSample}{"DRF"}=$h{DRF};
				$tmpF{$cSample}{"DRA"}=$h{DRA};
				$tmpF{$cSample}{"MQBP"}=$h{MQBP};
# 				$tmpF{$cSample}{"ICN"}=0;
				$tmpF{$cSample}{"FT"}=$h{FT};
								
				print "A:  $tmpF{$cSample}{freq} $h{DRA}\n" if $verbose;
			}
			
# 			next if average(\@fredS)<20;
			
			# assemble vcf line from CNVnator line # [$cSt,$cEn,$cLen,$cFreq,[$cSample],[$cSt],[$cEn],[$cLen],[$cFreq],[$cPval]]
			$CNVnatorID++;
			($cSt,$cEn,$cLen,$cFreq)=@$cL[0..3];
			$delOrDup=($cFreq<1)? "DEL":"DUP";
			print "B:  $cFreq $delOrDup $cChr:$cSt-$cEn\n"  if $verbose;
			
			$avgPVal=($avgPVal/$avgPValC);
			@outArr=($cChr,$$cL[0],"CNVnator_$CNVnatorID","N","<$delOrDup>",fredQ($avgPVal),"."); # ("CHROM","POS","ID","REF","ALT","QUAL","FILTER"); 
			
			# exclude CNR=$cFreq;
			push @outArr, "SVTYPE=$delOrDup;SVLEN=$cLen;END=$cEn;CHR2=$cChr;IMPRECISE;TOOL=CNVnator;EVTYPE=RD;LOCATION=$cChr:$cSt-$cEn;CNV=1;MQ=".round(average(\@MQS),1).";PCSD=".round(average(\@PCSD),1).";CPE=".round(average(\@CPE),0).";CSR=".round(average(\@CSR),0);
			push @outArr, join(":",@outFTypes);
			
			foreach $cSample (@$aSamples){
				
				$cGTS="";
				if(!exists($tmpF{$cSample})){
					$cGTS="0/0";
					for(1..(scalar(@outFTypes)-1)){ $cGTS.=":.";  }
				}else{
					$tmpF{$cSample}{"CNRD"}=$tmpF{$cSample}{"freq"};
					$tmpF{$cSample}{"CNP"}=fredQ($tmpF{$cSample}{"cPval"});
					
					$tmpF{$cSample}{"GT"}="1/1";
					#adjust genotype to het 0/1 if CNV
				  	if($delOrDup eq "DEL" and $tmpF{$cSample}{DRA} >=0.2 and $tmpF{$cSample}{DRA} <0.8){
						$tmpF{$cSample}{"GT"}="0/1"; 
					}elsif($delOrDup eq "DUP" and $tmpF{$cSample}{DRA} >1.2 and $tmpF{$cSample}{DRA} <=1.75){
						$tmpF{$cSample}{"GT"}="0/1"; 
					}
					
					map {  $cGTS.=$tmpF{$cSample}{$_}.":"  } @outFTypes; chop($cGTS);
				}
				push @outArr, $cGTS;
			}
			print OUT join("\t",@outArr)."\n";
			
			
		}
	}
}



sub ovlCNVnatorLine{

	($cChr,$cSt,$cEn,$SVTYPE,$cLen)=@_;
	return (-1) if !exists($iCNVmaster{$cChr});
	
	($iSt,$iEn)=binsearch_range(\@{$iCNVmaster{$cChr}},($cSt-$cLen),$cEn);
	return (-1) if $iSt<0; # no overlap
	
	$aChr=$CNVmaster{$cChr}; # pointer to current chromosome array 
# 	print "$cSt,$cEn : $iSt,$iEn : $cChr,$cSt,$cEn,$SVTYPE,$cLen\n";
	for my $i ($iSt..$iEn){
		next if exists($CNVMNatorLineBlackList{$cChr}{$i}); # already merged to lumpy
		$cL=$$aChr[$i]; # pointer to current line
		($cSt2,$cEn2,$cLen2,$cFreq2)=@$cL[0..3];
		$cOvlLen=overlap($cSt,$cEn,$cSt2,$cEn2);
		
		$lumpOvlDiff=abs($cLen-$cOvlLen);
		$cnvnOvlDiff=abs($cLen2-$cOvlLen);
		$lumpOvlPc=sprintf("%.1f",($cOvlLen/$cLen));
		$cnvnOvlPc=sprintf("%.1f",($cOvlLen/$cLen2)); # % of DGV variant
		
# 		print "ovl: l:$lumpOvlDiff c:$cnvnOvlDiff , l:$lumpOvlPc c:$cnvnOvlPc, abs Len: l:$cLen, c:$cLen2 \n";
		
# 		print "ovl:$cOvlLen $cSt,$cEn,$cLen,$SVTYPE  CNVnator $cSt2,$cEn2,$cLen2,$cFreq2 \n";
		if( ($lumpOvlPc >= $ovlCutOffLumpy or $lumpOvlDiff<=$ovlToleranceLumpy)
		and ($cnvnOvlPc >= $ovlCutOffLumpy or $cnvnOvlDiff<=$ovlTolerance)  
		and (($cFreq2>1.2 and $SVTYPE eq "DUP") or ($cFreq2<0.8 and $SVTYPE eq "DEL"))    ){
# 		if( $lumpOvlPc >= $ovlCutOffLumpy  and ($cOvlLen/$cLen2) >= $ovlCutOffLumpy ){
# 			print  "mergeToCluster $i:  $lumpOvlPc >= $ovlCutOffLumpy  and ($cOvlLen/$cLen2) >= $ovlCutOffLumpy and (($cFreq2>1.2 and $SVTYPE eq DUP) or ($cFreq2<0.8 and $SVTYPE eq DEL)) \n";
			return $i;
		}
	}
	return (-1);
}

sub binsearch_range{ # examples 201,600->300-600; 200,200->200; 150,160->-1; 40-99->-1, 1001,1001->-1; 50,100->100, 1000,1000->1000, 1000,1001->1000, get all st and with same coord
 my ($r,$sSt,$sEn)=@_; ($sSt,$sEn)=sort {$a <=> $b} ($sSt,$sEn);
 if(($sSt<$$r[0] and $sEn<$$r[0]) or ($sSt>$$r[$#$r] and $sEn>$$r[$#$r])){ return ((((-1),(-1)))); } # St and En outside array
 $iSt = binsearch_pos {$a <=> $b} $sSt, @$r;
 $iEn = binsearch_pos {$a <=> $b} $sEn, @$r;
 $iEn=$#$r if $iEn>$#$r; $iEn-- if $sEn<$$r[$iEn] and $iEn>0;
 while(($iEn+1)<=$#$r and $$r[$iEn]==$$r[$iEn+1]){ $iEn++; } # correct for multiple same entries upper end
 while(($iSt-1)>=0 and $$r[$iSt]==$$r[$iSt-1]){ $iSt--; } # correct for multiple same entries lower end
 if($iSt>$iEn){ return ((((-1),(-1))));} # 140,150 no interval in between
 return ($iSt,$iEn);
}


sub average{ my ($a)=@_; return 0 if @$a==0; $avgS=0; $avgC=0; foreach (@$a){ $avgC++; $avgS+=$_; } return round($avgS/$avgC,2)}
sub overlap{ my($R1,$R2,$Q1,$Q2)=@_; ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;}


sub log10 {
        my $n = shift;
        return log($n)/log(10);
}

sub fredQ{
	($pVal)=@_;
	$cQualS=($pVal==0)? 40:sprintf("%.0f",((-10)*log10($pVal)));
	$cQualS=40 if $cQualS>40;
	$cQualS=1 if $cQualS<=0;
	return $cQualS;
}

sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}

