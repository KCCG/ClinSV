
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use POSIX;
use Bio::DB::Big;
use List::Util qw(sum);
use File::Basename;

##### v4 
# - use bigwig bins, and 
# - add read depth ratio for MQ=60 portion, if portion > 30%
# - convert DEL or dups not passing the depth filter to BND, do not exclude

$readLen=150;

$inVCF=shift(@ARGV);
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

	# default to grch37 naming convention
	$chrPf='';
	$Y_chr="Y";
	$X_chr="X";
	$M_chr="MT";
}

%tabixC;
$tabixC{SR} = new Tabix(-data => "$jobfs/SR.brkp.gz",-index => "$jobfs/SR.brkp.gz.tbi");
$tabixC{PE} = new Tabix(-data => "$jobfs/PE.brkp.gz",-index => "$jobfs/PE.brkp.gz.tbi");

$typeCat{"DEL"}={ "1mD,2mU",1, "1pU,2pD",1, "mD,pU",1 };
$typeCat{"DUP"}={ "1mU,2mD",1, "1pD,2pU",1, "mU,pD",1 };
$typeCat{"INV"}={ "1mU,2pD",1, "1pD,2mU",1, "1pU,2mD",1, "1mD,2pU",1, "pD,pU",1, "mD,mU",1 };


$verbose=0;

if($inVCF=~/^(.+)\.vcf(.gz)*$/){ $outVCF=$1 } $outVCF.=".out";

print STDERR "# running add-depth-to-PE-SR-calls-v4.pl \n";
print STDERR "# in VCF: $inVCF\n";
print STDERR "# out VCF: $outVCF\n";
print STDERR "# in SR control: $jobfs/SR.brkp.gz\n";
print STDERR "# in PE control: $jobfs/PE.brkp.gz\n";


print STDERR "create tabix vcf to extract the copy neutral deletion duplications...\n";
############################################################


if(! -e "$inVCF.gz.tbi"){
`sort_bgzip $inVCF; tabix -f -p vcf $inVCF.gz`;
}

my $tabix = new Tabix(-data => "$inVCF.gz",-index => "$inVCF.gz.tbi");


print STDERR "open in the input VCF...\n";
############################################################


my $V=My::VCF2->new($inVCF);
$W=My::VCF2->new($inVCF);
$X=My::VCF2->new($inVCF);
foreach $cH (($V,$W,$X)){
$cH->add_header("##FORMAT=<ID=DRF,Number=1,Type=Float,Description=\"Read depth ratio of variant vs flanking regions\">");
$cH->add_header("##FORMAT=<ID=DRA,Number=1,Type=Float,Description=\"Read depth ratio of variant vs the average genome wide coverage\">");
$cH->add_header("##FORMAT=<ID=DBP,Number=1,Type=Integer,Description=\"Read depth at break point\">");
$cH->add_header("##FORMAT=<ID=MQBP,Number=1,Type=Integer,Description=\"Average read mapping quality of reads supporting both breakpoints\">");
$cH->add_header("##FORMAT=<ID=ICN,Number=1,Type=Integer,Description=\"Is the structural variant not a CNV? 1 = not, 0 = CNV. Yes if DRA or DRF smaller 0.8 or greater 1.2\">");
$cH->add_header("##FORMAT=<ID=IDD,Number=1,Type=Integer,Description=\"Is deletion compensated by duplication or vice versa: 1 = yes, 0 = no\">");
$cH->add_header("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Automated call confidence filter column. Values LOW, PASS or HIGH\">");
$cH->add_header("##INFO=<ID=LOCATION,Number=.,Type=String,Description=\"Genomic location (chr:start-end)\">");
$cH->add_header("##INFO=<ID=CNV,Number=.,Type=Integer,Description=\"Is the structural variant a CNV? 1 = yes, 0 = no. Yes if DRA or DRF smaller 0.8 or greater 1.2\">");
$cH->add_header("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Average read mapping quality of the variant\">");
$cH->add_header("##INFO=<ID=PCSD,Number=1,Type=Float,Description=\"Population coverage standard deviation of contorl cohort\">");

$cH->add_header("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome name of second breakpoint\">");
$cH->add_header("##INFO=<ID=CPE,Number=1,Type=Float,Description=\"Number discordant pairs in control cohort concordant to this variant \">");
$cH->add_header("##INFO=<ID=CSR,Number=1,Type=Float,Description=\"Number split reads in control cohort concordant to this variant \">");
$cH->add_header("##INFO=<ID=OL,Number=1,Type=Float,Description=\"Number links with at least 10 lines evidence to other genomic location using split reads and discordant pairs from control set.\">");
$cH->add_header("##FORMAT=<ID=PAFSU,Number=1,Type=Float,Description=\"Population variant allele frequency estimated from sum of discordant pairs (DP) and split reads (SR) in control cohort. \">");


}

open(VCF, ">$outVCF") || die "$outVCF nicht gefunden";
print VCF ${$V->print_header};

@aSamples=@{$$V{samples}};
print STDERR "samples in vcf: ".join(",",@aSamples)."\n";


print STDERR "open the bw files...\n";
############################################################
our %bwObj;

foreach $cSample (@aSamples){ 
	
	foreach $cT (("mq","q0","q20")){
				
# 		$cQBW=$jobfs."/".$cSample.".$cT.bw";
		$cQBW=$projectDir."/alignments/$cSample/bw/".$cSample.".$cT.bw";
		die " $cQBW ! exists" if (! -f $cQBW);
		$bwObj{$cSample}{$cT} = Bio::DB::Big->open($cQBW);
	}	
	print STDERR "# sample: $cSample, bw found\n";
}

$bwObj{mq} = Bio::DB::Big->open($S_control_bw_folder."/aMQ.bw");
$bwObj{sd} = Bio::DB::Big->open($S_control_bw_folder."/popCovStdev.bw");

#read in chr2len
open(IN1, "<$inRef.fai") || die "can not open $inRef";
while(<IN1>){ chomp; @_=split("\t",$_); $chr2len{$_[0]}=$_[1]; }close(IN1);


print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
############################################################
%meanCov;
foreach $cSample (@aSamples){ 
	$sumCov=0;
	for(1..22){
		$sumCov+=getCov($chrPf.$_,20000001,30000000,"q0",58,$cSample);
	}
	foreach $cChr ((keys %chr2len,"A")){ $meanCov{$cChr}{$cSample}=sprintf("%.2f",($sumCov)/22); }
	
	# determine X and Y coverage and round to 0, 0.5, 1, 1.5, 2, ... times average coverage
	$meanCov{$X_chr}{$cSample}=sprintf("%.0f",(getCov($X_chr,1,$chr2len{$X_chr},"q0",58,$cSample)/$meanCov{A}{$cSample})/0.5)*0.5*$meanCov{A}{$cSample};
	$meanCov{$Y_chr}{$cSample}=sprintf("%.5f", (getCov($Y_chr,6641419,10079253,"q0",0,$cSample)+getCov($Y_chr,13800704,23668908,"q0",0,$cSample))/2 ); # two mq 60 region one of each Y arm 
	print STDERR "Y:".$meanCov{$Y_chr}{$cSample}."\n";
	$meanCov{$Y_chr}{$cSample}=0 if $meanCov{$Y_chr}{$cSample} == (-1); # if -1 then because it is not covered	
	$meanCov{$Y_chr}{$cSample}=sprintf("%.0f",($meanCov{$Y_chr}{$cSample}/$meanCov{A}{$cSample})/0.5)*0.5*$meanCov{A}{$cSample};
	
	# MT coverage
	$meanCov{$M_chr}{$cSample}=sprintf("%.2f",(getCov($M_chr,1,$chr2len{$M_chr},"q0",58,$cSample)));	

	print STDERR "average coverage $cSample autosome:".$meanCov{A}{$cSample}.", X:".$meanCov{$X_chr}{$cSample}.", Y:".$meanCov{$Y_chr}{$cSample}.", MT:".$meanCov{$M_chr}{$cSample}."\n";
}


print STDERR "parse VCF...\n";
############################################################
############################################################
############################################################

$variantsAffectedByGene=0;
$allVariantsInVCF=0;

while ($V=$V->next_line){
  		
  	
	$allVariantsInVCF++;
  	$$V{INFO}{SR}=0 if !exists($$V{INFO}{SR});		
	$$V{INFO}{PE}=0 if !exists($$V{INFO}{PE});
	$$V{INFO}{SVLEN}=abs($$V{INFO}{SVLEN}) if exists($$V{INFO}{SVLEN});
	
	($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{INFO}{"END"});
	

	############################ create SVLEN for BNDs on same Chr
	
	if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			($cChr2,$cPos2)=($1,$2);
			$$V{INFO}{"END"}=$cPos2;
			
			$$V{INFO}{CHR2}=$cChr2;
			if($$V{CHROM} eq $cChr2){
				$$V{INFO}{SVLEN}=abs($$V{POS}-$cPos2);
				($p1,$p2) = sort {$a <=> $b} ($$V{POS}, $cPos2); 
				$$V{INFO}{LOCATION}=$$V{CHROM}.":".$p1."-".$p2;
			}else{
				$$V{INFO}{LOCATION}=$$V{CHROM}.":".$$V{POS}.",".$cChr2.":".$cPos2;

			}		
			
	}else{
		$$V{INFO}{LOCATION}=$$V{CHROM}.":".$$V{POS}."-".$$V{INFO}{"END"};
		($cChr2,$cPos2)=($$V{CHROM},$$V{INFO}{"END"});
		$$V{INFO}{CHR2}=$$V{CHROM};
	}
	
	next if $$V{CHROM} eq "hs37d5" or $$V{INFO}{CHR2} eq "hs37d5";
	


# 	next if $$V{INFO}{SVTYPE} ne "BND";
	
	next if uc($cChr) eq $Y_chr and $meanCov{$Y_chr}{$cSample}==0 and $$V{INFO}{SVTYPE} ne "DUP" and $$V{INFO}{SVTYPE} ne "BND"; # ignore DEL,INV,BND (not DUP) on Y chromosome if Y has coverage == 0;

# 	print "### in: ".${$V->print_line};
# 	print "1 \n";	
	############################ pre filter remove low evidence SVs
	$cPASS_C=0;
	foreach $cSample (@aSamples){
		$cWeightA=($$V{GTS}{PE}{$cSample}+$$V{GTS}{SR}{$cSample}*2);
		$cWeightB=($$V{GTS}{PE}{$cSample}+$$V{GTS}{SR}{$cSample});
		if ($cWeightA>=3){
			
			($cDRF,$cDRA)=(-1,-1);
			calcDRF($V,$cSample);
			$$V{GTS}{DRF}{$cSample}=$cDRF;
			$$V{GTS}{DRA}{$cSample}=$cDRA;
			
			
			$$V{GTS}{ICN}{$cSample}=(($cDRF>=0.8 and $cDRF<=1.2) or ($cDRA>=0.8 and $cDRA<=1.2))? 1:0;
			if ($cWeightB>=4 or ($cDRF<0.8 and $$V{INFO}{SVTYPE} eq "DEL") or ($cDRF>1.2 and $$V{INFO}{SVTYPE} eq "DUP")){
				$cPASS_C++;
			}
		}
	}
# 	print "2 \n";
# 	print ${$V->print_line} if $cPASS_C==0;

# 	print "### skip: not passing pre-filter ".${$V->print_line} if $verbose and $cPASS_C==0;
	next if $cPASS_C==0;
	
# 	print "3 \n";
	############################ Assign PASS filter
	map { assignPASS($V,$_) } @aSamples;
	undef @pSamples; map {push @pSamples,$_ if $$V{GTS}{FT}{$_} eq "PASS"} @aSamples;
	print "### skip: not passing ".${$V->print_line} if $verbose and @pSamples==0;
# 	print ${$V->print_line} if @pSamples==0;
	next if @pSamples==0;
	
	############################ check if deldupe evidence exists
  	# if event is overlapping check for coverage dip at non intersecting parts
# 	print "4 \n";
  	checkIfDelDup($V) if $$V{INFO}{SVTYPE} eq "DEL" or $$V{INFO}{SVTYPE} eq "DUP";	
# 	print "5 \n";


	
	############################ output samples with copy neutral deletions or duplications seperatly, add info flag ICN 

	# how many concordant pairs are there per SV in the background dataset, exclude CNVs
	# how many other predominant destinations are there (discard cluster <10)
# 	if($$V{INFO}{SVTYPE} eq "BND" or $$V{INFO}{SVTYPE} eq "INV" or $$V{INFO}{CNV}==1){	
	concordantSRPEs($V,$$V{CHROM},$$V{POS},$cChr2,$cPos2,$$V{INFO}{SVTYPE}); # adds ($$V{INFO}{CDISC},$$V{INFO}{OLINK})

# 	print "6 \n";
	$cSUC=($$V{INFO}{CPE}+$$V{INFO}{CSR});
	
	foreach $cSample (@aSamples){
		$cSU=$$V{GTS}{SR}{$cSample}+$$V{GTS}{PE}{$cSample};
				
		$$V{GTS}{PAFSU}{$cSample}=($cSU>0)? roundA(  ($cSUC/$cSU)/$S_SV_control_numSamples   ):(-1);  
# 			print " ${$breaks{$cBrID}}[3] $cSUC / $cSU ) / $S_SV_control_numSamples = ".${$BrStat{PAFSU}}[$cID]."\n";
	} 
	

	######## split by type: CN DUP, CN DEL, copy neutral (as is)
	# SVTYPE BND is kept 
	if($$V{CHROM} eq $cChr2){
		undef %sample2CN; undef %CN2Sample;
	
		foreach $cSample (@aSamples){
			
			if($$V{GTS}{ICN}{$cSample}==1){ $cCN="NEUTR" }
			elsif( $$V{GTS}{DRA}{$cSample}>1.2 ){ $cCN="DUP" }
			elsif( $$V{GTS}{DRA}{$cSample}<0.8 and $$V{GTS}{DRA}{$cSample}>=0 ){ $cCN="DEL" }
			elsif( $$V{GTS}{DRF}{$cSample}>1.2 ){ $cCN="DUP" }
			elsif( $$V{GTS}{DRF}{$cSample}<0.8 and $$V{GTS}{DRF}{$cSample}>=0 ){ $cCN="DEL" }
			else{ $cCN="."}
			push @{$CN2Sample{$cCN}}, $cSample if $$V{GTS}{FT}{$cSample} ne "FAIL";
			$sample2CN{$cSample}=$cCN;
		}
	
		$idOrig=$$V{ID};
		
		if ($idOrig =~ /_[12]$/){ chop($idOrig); }
		else{$idOrig.="_";}
# 		print " $idOrig $$V{ID} \n";

		for $cGorL (("DUP","DEL")){ # output DUP and DELs seperatly
			if( @{$CN2Sample{$cGorL}} > 0 ){ # if at least 1 is CN DEL/DUP than output it as seperate line
				$W=$W->parse_line($V->print_line);
				$$W{INFO}{CNV}=1;
				$$W{INFO}{SVTYPE}=$cGorL;
				foreach $cSample (@aSamples){
					if($sample2CN{$cSample} eq $cGorL){ # mask DUP/DEL samples in orig $V (
						$$V{GTS}{FT}{$cSample}="FAIL";
						#foreach $cGTS (keys %{$$V{GTS}}){ $$V{GTS}{$cGTS}{$cSample}="." if exists $$V{GTS}{$cGTS}{$cSample}; }
					}else{ # mask NEUTR samples in new $W (the CNV)
						$$W{GTS}{FT}{$cSample}="FAIL";
						#foreach $cGTS (keys %{$$W{GTS}}){ $$W{GTS}{$cGTS}{$cSample}="." if exists $$W{GTS}{$cGTS}{$cSample}; }
					}
				}
				
				
				# in case lumpy called a DEL but the CN is gain that make it a BND
				$cGainOrLoss=($cGorL eq "DEL")? "loss":"gain"; 
				$id_stem=$idOrig.$cGainOrLoss;
				
				next if exists($$V{INFO}{SECONDARY}); # for BNDs skip the secondary, and write it by copying primary
				
				if($cGorL ne $$V{INFO}{SVTYPE} ){
					
					($strand,@aStrands)=getStrand($V);
# 					print "sStrands: $sStrands @aStrands\n";
					if( $strand eq '+-' ){
						
						$$W{ID}=$id_stem."_1";
						$$W{INFO}{SVTYPE}="BND";  
						$$W{ALT}="N[".$$W{CHROM}.":".$$W{INFO}{"END"}."[";
						$$W{INFO}{STRANDS}="+-:".$strand_arr[1];  
						$$W{INFO}{MATEID}=$id_stem."_2";
						print VCF ${$W->print_line};
						
						$$W{ID}=$id_stem."_2";
						$$W{ALT}="]".$$W{CHROM}.":".$$W{POS}."]N";
						$$W{POS}=$$W{INFO}{"END"};
						$$W{INFO}{STRANDS}="-+:".$strand_arr[1];  
						$$W{INFO}{MATEID}=$id_stem."_1";
						$$W{INFO}{SECONDARY}=1;
						print VCF ${$W->print_line};
						
					
					}elsif( $strand eq '-+' ){
						
						$$W{ID}=$id_stem."_1";
						$$W{INFO}{SVTYPE}="BND";  
						$$W{ALT}="]".$$W{CHROM}.":".$$W{INFO}{"END"}."]N";
						$$W{INFO}{STRANDS}="-+:".$strand_arr[1];
						$$W{INFO}{MATEID}=$id_stem."_2";
						print VCF ${$W->print_line};
						
						$$W{ID}=$id_stem."_2";

						$$W{ALT}="N[".$$W{CHROM}.":".$$W{POS}."[";
						$$W{POS}=$$W{INFO}{"END"};
						$$W{INFO}{STRANDS}="+-:".$strand_arr[1];  
						$$W{INFO}{MATEID}=$id_stem."_1";
						$$W{INFO}{SECONDARY}=1;
						print VCF ${$W->print_line};
						
					}elsif( $strand eq '++' ){
						
						$$W{ID}=$id_stem."_1";
						$$W{INFO}{SVTYPE}="BND";  
						$$W{ALT}="N]".$$W{CHROM}.":".$$W{INFO}{"END"}."]";
 						$$W{INFO}{MATEID}=$id_stem."_2";
						print VCF ${$W->print_line};
						
						$$W{ID}=$id_stem."_2";
						$$W{ALT}="N]".$$W{CHROM}.":".$$W{POS}."]";
						$$W{POS}=$$W{INFO}{"END"};

						$$W{INFO}{MATEID}=$id_stem."_1";
						$$W{INFO}{SECONDARY}=1;
						print VCF ${$W->print_line};
			
					}elsif( $strand eq '--' ){
						
						$$W{ID}=$id_stem."_1";
						$$W{INFO}{SVTYPE}="BND";  
						$$W{ALT}="[".$$W{CHROM}.":".$$W{INFO}{"END"}."[N";
 						$$W{INFO}{MATEID}=$id_stem."_2";
						print VCF ${$W->print_line};
						
						$$W{ID}=$id_stem."_2";
						$$W{ALT}="[".$$W{CHROM}.":".$$W{POS}."[N";
						$$W{POS}=$$W{INFO}{"END"};

						$$W{INFO}{MATEID}=$id_stem."_1";
						$$W{INFO}{SECONDARY}=1;
						print VCF ${$W->print_line};
						
					}else{
						die "unexpected STRAND: \"$strand\" ".${$W->print_line}."\n";
					}
				
				}else{
					$$W{ID}=$id_stem;
					print "### out $cGorL: ".${$W->print_line} if $verbose;
					print VCF ${$W->print_line};
				}
				

# 				print "### out $cGorL: ".${$W->print_line} if $verbose;
				
				
				
			}
		}
	}
	
	
	
	
	$$V{INFO}{CNV}=0;
	print "### out NEUTR: ".${$V->print_line} if $verbose;
	print VCF ${$V->print_line} if scalar(grep { $$V{GTS}{FT}{$_} eq "."  } @aSamples) != @aSamples and scalar(grep { $$V{GTS}{FT}{$_} eq "FAIL"  } @aSamples) != @aSamples;

			
		
	
# 	}

	
}


sub assignPASS {

	my ($V,$cSample)=@_;
	($cChr,$cSt,$cChr2,$cEn)=($$V{CHROM},$$V{POS},$$V{INFO}{CHR2},$$V{INFO}{"END"});
	$cLen=($cEn-$cSt+1);
	$cType=$$V{INFO}{SVTYPE};
	
# 	print STDERR "1.1 ";
	
	if (!exists($$V{GTS}{DRF}{$cSample}) or !exists($$V{GTS}{DRA}{$cSample})){
		($cDRF,$cDRA)=(-1,-1);
		calcDRF($V,$cSample) ;
		$$V{GTS}{DRF}{$cSample}=$cDRF;
		$$V{GTS}{DRA}{$cSample}=$cDRA;
	}
	
	$cDRF=$$V{GTS}{DRF}{$cSample}; # depth ratio flanking
	$cDRA=$$V{GTS}{DRA}{$cSample}; # depth ratio to average
	
	$$V{GTS}{ICN}{$cSample}=(($cDRF>=0.8 and $cDRF<=1.2) and ($cDRA>=0.8 and $cDRA<=1.2))? 1:0;
	
	$$V{GTS}{DBP}{$cSample}=calcDBP($V,$cSample); # depth at breakpoint
# 	print STDERR "1.2 ";
	$$V{GTS}{MQBP}{$cSample}=calcMQBP($V,$cSample); # mapping quality at breakpoint
# 	print STDERR "1.3 ";
	
	if($cChr2 eq $cChr){
		$$V{INFO}{MQ}=getMQSD($cChr,$cSt,$cChr2,$cEn,$cSample,"mq") if !exists($$V{INFO}{MQ});
# 		print STDERR "1.4 ";
		$$V{INFO}{PCSD}=getMQSD($cChr,$cSt,$cChr2,$cEn,$cSample,"sd") if !exists($$V{INFO}{PCSD});
# 		print STDERR "1.4.1 ";
		if( $$V{INFO}{MQ} == (-1) or $$V{GTS}{DRA} <0.2 ){  # if unsuccessfull or homozygous deletion use flanking regions to determine MQ
			$cMQ1=getMQSD($cChr,($cSt-$cLen),$cChr,($cSt-1),$cSample,"mq");
			$cMQ2=getMQSD($cChr2,($cEn+1),$cChr2,($cEn+$cLen),$cSample,"mq");
			$$V{INFO}{MQ}=round(avgIfGr0($cMQ1,$cMQ2),1);	
		}
# 		print STDERR "1.5 ";
	}
	
# 	print STDERR "1.6 \n";
	
	$cICN=$$V{GTS}{ICN}{$cSample}; # is copy neutral
	$cDBP=$$V{GTS}{DBP}{$cSample}; # depth at breakpoint
	$cPE=$$V{GTS}{PE}{$cSample};
	$cSR=$$V{GTS}{SR}{$cSample};
	
	
	
	#### trust all events but deletions <1000 bases if at least 4 lines of evidence are present
	# Deletions <1000 can be an artifacts from local high density of discordant PEs. They definitely need more evidence
	# duplication have a different orientation so there are automatically less false positives
	if ($cType eq "BND" or $cType eq "INV" or $cType eq "DUP" 
	or ($cType eq "DEL" and $$V{INFO}{SVLEN}>=1000)   ){
	
		if(($cPE+$cSR)>=4){ $$V{GTS}{FT}{$cSample}="PASS"; return }
		else{ $$V{GTS}{FT}{$cSample}="FAIL"; return }
	}
		
		
	### SR and PE alone
	$$V{GTS}{FT}{$cSample}="PASS",return if( $cSR>=2 & $cPE>=2 & ($cSR+$cPE)>=6 & $cLen<=10000 & $cDBP<150 );
	$$V{GTS}{FT}{$cSample}="FAIL",return if ( $cSR+$cPE )==0;
	
	
	### the following criteria need depth ratio, if it could not be determined, then FT=FAIL
	if($cDRF < 0 and $cDRA < 0){  $$V{GTS}{FT}{$cSample}="FAIL"; return  }
	
	$cDSC=0; # if depth could not be 
	if($cType eq "DEL"){ if($cDRF<0.8 or $cDRA<0.8){$cDSC+=2} if($cDRF<0.2 or $cDRA<0.2){$cDSC+=2} }
	
	undef %cScores;
	$cScores{SVLEN}=sprintf("%.1f",(($cLen/84)/1.6+1)); # 10=1.1, 100=1.7, 200=2.5, 500=4.7, 800=7.0
	$cScores{PE}=($cPE/2+3); # PE 1=3.5, 2=4, 3=4.5, 4=5, 5=5.5, 6=6, 7=6.5, 8=7, 9=7.5, 10=8
	$cScores{SR}=($cSR/2+4); # PE 1=4.5, 2=5, 3=5.5, 4=6, 5=6.5, 6=7, 7=7.5, 8=8, 9=8.5, 10=9
	if($cType eq "DEL"){ $cScores{"cov"}=($cDRF*(-10)+13); $cScores{"cov"}=($cDRA*(-10)+13) if ($cDRA*(-10)+13)>$cScores{"cov"}; } # DRF 0.8=5, 0.7=6, 0.6=7, 0.5=8, 0.4=9, 0.3=10, 0.2=11, 0.1=12, 0.1=13
	
	foreach $cScT ( keys %cScores ){ # each score is allowed to be 10 at most
		$cScores{$cScT}=($cScores{$cScT}>10)? 10:$cScores{$cScT};
		$cScores{$cScT}=($cScores{$cScT}<0)?  0:$cScores{$cScT};
	}
	
# 	$$V{INFO}{CIEND95}=$cScores{SVLEN}.",".$cScores{PE}.",".$cScores{SR}.",".$cScores{"cov"};
# 	$$V{INFO}{CIEND}=($cScores{SVLEN}+$cScores{PE}+$cScores{SR}+$cScores{"cov"});
# 	$$V{GTS}{FT}{$cSample}="PASS"; return;
	
	### SR and depth alone
	$$V{GTS}{FT}{$cSample}="PASS",return if $cSR>=2 & $cDSC>=2 & $cDBP<150 & ($cScores{SVLEN}+$cScores{SR}+$cScores{"cov"}) > 17;
	
	### PE and depth alone	
	$$V{GTS}{FT}{$cSample}="PASS",return if $cPE>=3 & $cDSC>=2 & $cDBP<150 & ($cScores{SVLEN}+$cScores{PE}+$cScores{"cov"}) > 17;
	
	### SR and PE and depth
	$$V{GTS}{FT}{$cSample}="PASS",return if( $cSR>=1 & $cPE>=1 & ($cSR+$cPE)>=4 & $cDSC>=2 & $cDBP<150);
	$$V{GTS}{FT}{$cSample}="PASS",return if( $cSR>=2 & $cPE>=2 & ($cSR+$cPE)>=8 & $cDSC>=2);
	
	$$V{GTS}{FT}{$cSample}="FAIL";

}


sub getStrand {

	($V)=@_;
	
	@strand_tmp=split(":",$$V{INFO}{STRANDS});
	$sStrands=$strand_tmp[0];
	@aStrands=split("",$strand_tmp[0]);
	
	
	return ($sStrands,@aStrands);
}


sub getBps {

	($V)=@_;
	
	($sStrands,@aStrands)=getStrand($V);
	
	($sChr,$sPos,$sOri,$eChr,$ePos,$eOri)=($$V{CHROM},$$V{POS},$aStrands[0],$$V{INFO}{CHR2},$$V{INFO}{"END"},$aStrands[1]);
	
	if ($sChr eq $eChr and $sPos>$ePos){
		@bps=({'chr',$eChr,'pos',$ePos,'ori',$eOri},{'chr',$sChr,'pos',$sPos,'ori',$sOri} );
	}else{
		@bps=({'chr',$sChr,'pos',$sPos,'ori',$sOri},{'chr',$eChr,'pos',$ePos,'ori',$eOri} );
	}
	
# 	print "  $$V{INFO}{SVTYPE}  'chr',$sChr,'pos',$sPos,'ori',$sOri   'chr',$eChr,'pos',$ePos,'ori',$eOri \n";
	
	
	return @bps;
}

sub calcMQBP {

	($V)=@_;
	@bps=getBps($V);
	
	@sum=();
	
	foreach (@bps){
		if ($$_{ori} eq "+"){
# 			print $$_{chr}." ".($$_{pos}-$readLen)." ".$$_{chr}." ".($$_{pos}-1)." $readLen \n";
			push @sum, getMQSD($$_{chr},($$_{pos}-$readLen),$$_{chr},($$_{pos}-1),$cSample,"mq");
		}elsif ($$_{ori} eq "-"){
# 			print $$_{chr}." ".($$_{pos}+1)." ".$$_{chr}." ".($$_{pos}+$readLen)." $readLen\n";
			push @sum, getMQSD($$_{chr},($$_{pos}+1),$$_{chr},($$_{pos}+$readLen),$cSample,"mq");
		}
	}
	
# 	print " getMQSD :".join("\t",@sum)."\n";
# 	$tmp=<STDIN>;
	
	return round(sum(@sum)/scalar(@sum),0);	
# 	return round(avgIfGr0Arr(\@sum),0);	

}

sub calcDBP {


	($V,$cSample)=@_;
	@bps=getBps($V);
	
	@sum=();
	
	foreach (@bps){
		if ($$_{ori} eq "+"){
			push @sum, getCov($$_{chr},($$_{pos}-$readLen),($$_{pos}-1),"q0",0,$cSample);
		}elsif ($$_{ori} eq "-"){
			push @sum, getCov($$_{chr},($$_{pos}+1),($$_{pos}+$readLen),"q0",0,$cSample);
		}
	}
	
# 	print " getMQSD $cSample :".join("\t",@sum)."\n";
# 	$tmp=<STDIN>;
	
	return round(sum(@sum)/scalar(@sum),0);	
# 	return round(avgIfGr0Arr(\@sum)/2,0);	


}
		
sub calcDRF {

	($V,$cSample)=@_;
	
	@bps=getBps($V);
	
	if($bps[0]->{chr} ne $bps[1]->{chr}){
		$cDRF=(-1);
		$cDRA=(-1);	
		return
	}
	
	($cChr,$cSt,$cEn)=($bps[0]->{chr},$bps[0]->{pos}, $bps[1]->{pos});
	$cLen=($cEn-$cSt+1);
	
# 	print $bps[0]->{chr}."\n";

	
	die " cSt>cEn  $cSt>$cEn  "  if $cSt>$cEn;
	
	# for flanking only rely on high MQ, else take global mean
	$cDFL1=getCov($cChr,($cSt-$cLen),($cSt-1),"q20",50,$cSample);
	$cDFL2=getCov($cChr,($cEn+1),($cEn+$cLen),"q20",50,$cSample);
	$cDFL=avgIfGr0($cDFL1,$cDFL2);
	$cDFL=$meanCov{$cChr}{$cSample} if $cDFL<=0;
	
	# for internal try high MQ else all
	$cDIN=getCov($cChr,$cSt,$cEn,"q20",50,$cSample);
	$cDIN=getCov($cChr,$cSt,$cEn,"q20",0,$cSample) if $cDIN<0;
	print STDERR "$cSample: cDIN: $cDIN $cChr,$cSt,$cEn meanCov:$meanCov{$cChr}{$cSample} \n" if $verbose;
	if($cDIN<0){return}
	$cDRF=($cDFL<=0)? (-1):sprintf("%.2f",$cDIN/$cDFL);
	$cDRA=($meanCov{$cChr}{$cSample}<=0)? (-1):sprintf("%.2f",$cDIN/$meanCov{$cChr}{$cSample});
	
	$cDRF=(-1) if $cDRF<0;
	$cDRA=(-1) if $cDRA<0;

	return;
	
}




sub avgIfGr0Arr{
	my ($arr)=@_;
	
	$retVal=(-1);
	if($$arr[0]>0 and $$arr[1]>0){ $retVal=($$arr[0]+$$arr[1])/2;  }
	elsif($$arr[0]>0){ $retVal=$$arr[0];  }
	elsif($$arr[1]>0){ $retVal=$$arr[1];  }
	return $retVal;
}

sub avgIfGr0{
	my ($v1,$v2)=@_;
	$retVal=(-1);
	if($v1>0 and $v2>0){ $retVal=($v1+$v2)/2;  }
	elsif($v1>0){ $retVal=$v1;  }
	elsif($v2>0){ $retVal=$v2;  }
	return $retVal;
}
		
sub checkIfDelDup{
	
	# 1. check if evidence from counter call
	# 2. check evidence from coverage
	
	($V)=@_;
	undef @aOvl;
	die "!END before start: ".${$V->print_line} if $$V{INFO}{"END"}<$$V{POS};
	my $verbose=1;
	
	($qChr,$qPos1,$qPos2)=($$V{CHROM},$$V{POS},$$V{INFO}{"END"});
	$qLen=abs($qPos2-$qPos1);
	
	$ovlTolerance=($qLen<1000)? 0.5:0.25;
	$qiPos1=int($qPos1-($qLen*$ovlTolerance));
	$qiPos2=int($qPos1+($qLen*$ovlTolerance));
	return;
	$res = $tabix->query($qChr,$qiPos1,$qiPos2); 
	print "###??? checkCopyNeutral: ".$$V{INFO}{SVTYPE}." $qChr,$qPos1,$qPos2 len: $qLen, with tolerance: $qiPos1 $qiPos2\n" if $verbose;
# 	print STDERR "A $qChr,$qiPos1,$qiPos2\n";
	if(defined $res->get){
	while(my $line = $tabix->read($res)){
		
		$W=$W->parse_line(\$line);
		next if $$W{ID} eq $$V{ID};
		next unless ($$V{INFO}{SVTYPE} eq "DUP" and $$W{INFO}{SVTYPE} eq "DEL") or ($$V{INFO}{SVTYPE} eq "DEL" and $$W{INFO}{SVTYPE} eq "DUP");		
		
		($tChr,$tPos1,$tPos2)=($$W{CHROM},$$W{POS},$$W{INFO}{"END"});
		$tLen=abs($tPos2-$tPos1);
		
		$cOvl=overlap($qPos1,$qPos2,$tPos1,$tPos2);
		$qOvlFrac1=sprintf("%.5f",($cOvl/$qLen));
		$tOvlFrac2=sprintf("%.5f",($cOvl/$tLen));
		
		next unless ($qOvlFrac1>=0.5 and $tOvlFrac2>=0.5 and $qLen<1000) or ($qOvlFrac1>=0.75 and $tOvlFrac2>=0.75);
		
		foreach $cSample (@aSamples){ # the best del dup pair has to be selected separately for each sample
		
			$vPE=$$V{GTS}{PE}{$cSample};
			$vSR=$$V{GTS}{SR}{$cSample};
			$wPE=$$W{GTS}{PE}{$cSample};
			$wSR=$$W{GTS}{SR}{$cSample};
			next unless ($vSR>=1 and $vPE>=1) or ($vSR+$vPE)>=3;
			next unless ($wSR>=1 and $wPE>=1) or ($wSR+$wPE)>=3;
			$$V{GTS}{"IDD"}{$cSample}=1;
			
			print "-- IDD: $cSample $vPE:$vSR $wPE:$wSR\n" if $verbose;	
		}
		
		print "-- ovl: ".$$W{INFO}{SVTYPE}." $cOvl,$qOvlFrac1,$tOvlFrac2 \n" if $verbose;
		print "  tabix vcf lien: ".${$W->print_line} if $verbose;	
		
	}}
	
	return $V;
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
     

	my $aCov_stat = $bwObj{$cSample}{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	my $aMQ_stat = $bwObj{mq}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean') if $qMQ>0;	
	
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

sub getMQSD{
	my ($qChr,$qSt,$qChr2,$qEn,$cSample,$MQorSD)=@_;
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	$qSt=0 if $qSt<0;
	$qEn=$chr2len{$qChr} if $qEn > $chr2len{$qChr};
	return (-1) if $qChr ne $qChr2 or $qSt>=$qEn;
	
	$qSt=1 if $qSt<=0;
	$qEn=$chr2len{$qChr} if $chr2len{$qChr}<$qEn;
	
	my $qLen=($qEn-$qSt+1);
	my ($sumMQ,$fBins)=(0,0);
    my $qBins=100;

    $qBins=$qLen if ($qLen/$qBins)<1;   
    
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



sub overlap{
	my($R1,$R2,$Q1,$Q2)=@_; 
	($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	$ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	$returnVal=($ovlLen>0)? $ovlLen:0;
	return $returnVal;
}

sub median { ($a)=@_; @$a= sort { $a <=> $b } @$a; return $$a[((ceil(scalar(@$a)/2))-1)]; }



sub concordantSRPEs{


	my ($V,$cChr1,$cPos1,$cChr2,$cPos2,$cType)=@_;

	my @ol=();
	my $verbose=0;
	my $fl=1000;
	
	# adjust flanking if DEL <1k to encount better for random noise
	$fl=500 if ( $cType eq "DEL" and $cChr1 eq $cChr2 and abs($cPos1-$cPos2)<1000 );
	
	@rp = sort { $$a[0] cmp $$b[0] || $$a[1] <=> $$b[1];  } ([$cChr1,$cPos1],[$cChr2,$cPos2]);
	
	push @ol, [ $rp[0][0],$rp[0][1],$rp[1][0],$rp[1][1],1];
	
	$minSt=($rp[0][1]-$fl); $minSt=1 if $minSt<0;
	$maxSt=($rp[0][1]+$fl);	

	print STDERR "# $cChr1,$cPos1,$cChr2,$cPos2,$cType :: $rp[0][0],$minSt,$maxSt $$V{all}\n" if $verbose;
	
	
	foreach $cSRPE (("SR","PE")){

	$ci=0;
	$$V{INFO}{"C".$cSRPE}=0;
	
	
	$res = $tabixC{$cSRPE}->query($rp[0][0],$minSt,$maxSt);

	if(defined $res->get){$rC=0;
	while(my $line = $tabixC{$cSRPE}->read($res)){
		last if $rC++>$S_SV_control_numSamples*100;
		# 1       10100   type=1m,2p;type1=1m;NBP1=10167;MQ1=2;chr2=2;BP2=33141456;type2=2p;NBP2=33141533;MQ2=0	
		($bChr1,$bPos1,$bName)=split("\t",$line);
		undef %t; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $t{$cTMP2[0]}=$cTMP2[1]; }
		($bChr2,$bPos2,$bType)=($t{chr2},$t{BP2},$t{type});

				
		##### put into tmp hash for OL
		$clusteredYesNo=0;
		for($i=0; $i<=$#ol; $i++){
			($Bm1Chr,$Bm1Brkp,$Bm2Chr,$Bm2Brkp)=@{$ol[$i]};
			
			if (($bChr1 eq $Bm1Chr and $bChr2 eq $Bm2Chr) and # chr have to match
			   ($bPos1 < ($Bm1Brkp+$fl) and $bPos1 > ($Bm1Brkp-$fl)  ) and # pos 1 has to match
			   ($bPos2 < ($Bm2Brkp+$fl) and $bPos2 > ($Bm2Brkp-$fl)  )  ){ # pos 2 has to match	
			
				$ol[$i][4]++; $clusteredYesNo++; 			
# 				print STDERR "## $i $cSRPE ## clsuter $Bm1Chr,$Bm1Brkp,$Bm2Chr,$Bm2Brkp,  $bChr2=$Bm1Chr,  $bChr1=$Bm2Chr,  br1Diff:".(abs($Bm1Brkp-$bPos1)).", br2Diff:".(abs($Bm2Brkp-$bPos2))."  $ol[$i][4] \n" if $verbose;
				last;
			}
		}
		push @ol, [$bChr1,$bPos1,$bChr2,$bPos2] if $clusteredYesNo==0;
		## next if type is not satisfied
		
		##### next if location is not satisfied for CSE and PES
		next if $rp[1][0] ne $bChr2 or $bPos2<($rp[1][1]-$fl) or $bPos2>($rp[1][1]+$fl);
		
		
		if($cType eq "DEL"){		####### DEL, location test will result automatically in same locations
			#      74 type=1mD,2mU # DEL  <---   ---|
			#     511 type=1pU,2pD # DEL  |---   --->    
			#     299 type=mD,pU # DEL   --->   <--- 
# 			print STDERR "## hit $cSRPE ##$bChr1,$bPos1,$bChr2,$bPos2,$bType $cType \n" if $verbose;# and exists($typeCat{"DEL"}{$bType});

			$$V{INFO}{"C".$cSRPE}++ if exists($typeCat{"DEL"}{$bType});
			
		}elsif($cType eq "DUP"){    	####### DUP
			#     146 type=1mU,2mD # DUP   ---|   <---
			#       0 type=1pD,2pU # DUP   --->   |---
			#    3074 type=mU,pD # DUP   <---   ---> 
			$$V{INFO}{"C".$cSRPE}++ if exists($typeCat{"DUP"}{$bType});
		
		}elsif($cType eq "BND"){ 		######## BND can be any
			#    2243 type=1m,2m # ]p]N
			#    2393 type=1p,2m # N]p]    
			#     754 type=1m,2p # [p[N        
			#     721 type=1p,2p # N[p[
			$$V{INFO}{"C".$cSRPE}++;
			
		}elsif($cType eq "INV"){ 		####### INV
			#      88 type=1mU,2pD # INV   ---|   --->
			#      77 type=1pD,2mU # INV   <---   |---
			#     100 type=1pU,2mD # INV   |---   <---       
			#       0 type=1mD,2pU # INV   --->   ---|  
			#     194 type=pD,pU # INV   --->   --->
			#     146 type=mD,mU # INV   <---   <--- 
			$$V{INFO}{"C".$cSRPE}++ if exists($typeCat{"INV"}{$bType});

		}
		$ci++;
	}}

	print STDERR "# $cSRPE ".$$V{INFO}{"C".$cSRPE}."  \n" if $verbose;
	}
	
	# cluster per location, remove cluster with less than 10 lines of evidence
	# OL = number of remaining cluster 
	$$V{INFO}{OL}=scalar(@ol)-1;
	
	print STDERR "# OL ".$$V{INFO}{OL}." ".join("#",@{$ol[0]})."  ".join("#",@{$ol[1]})."  \n" if $verbose;
	
}





sub round{ my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec; }

sub roundA{ my ($cFreq)=@_;  if($cFreq!=1){ $RL=int(log10($cFreq)*(-1))+2; $cFreq=round($cFreq,$RL);  }  return $cFreq }
sub log10 { my $n = shift;  my $r=($n>0)? log($n)/log(10):0; return $r; }















