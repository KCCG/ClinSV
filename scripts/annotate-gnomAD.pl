
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::InterS qw(merge_segs);
use My::VCF2; 

use List::Util qw(min max);

($inBedGZ,$inVCF)=@ARGV;

%allowedChr;

$verbose=0;

############################# annotate the VCF file

BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;

my $tabix = new Tabix(-data => "$inBedGZ",-index => "$inBedGZ.tbi");


#definition flanking seg-dups: max 20% of len at +- 20% of variant ends



open(IN1, "zcat $inBedGZ |") || die " nicht gefunden";
while(<IN1>){ chomp; @_=split("\t",$_); 
	$allowedChr{$_[0]}++; 
	$c++; 
}close(IN1);

close(IN1);



my $VCFObj=My::VCF2->new($inVCF);

#$VCFObj->add_header("##INFO=<ID=PAC,Number=1,Type=Integer,Description=\"Population variant allele count in control cohort.\">");
$VCFObj->add_header("##INFO=<ID=PAFG,Number=1,Type=Integer,Description=\" Population variant frequency gnomAD.\">");

print ${$VCFObj->print_header};

while (my $V=$VCFObj->next_line){
  	
  	
  	$allVariantsInVCF++;

  	
  	($cChr1,$cPos1)=($$V{CHROM},$$V{POS});
	
	#$verbose=$$V{"ID"} eq "Lumpy_10987_loss" ? 1:0;
	
	if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
# 		next if exists($$V{"INFO"}{"SECONDARY"});
		if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			
			($cChr2,$cPos2)=($1,$2);
			if($cChr1 eq $cChr2){ 	
				$$V{"INFO"}{"SVLEN"}=abs($cPos1-$cPos2);
				$$V{"INFO"}{"END"}=$cPos2;
			}
		}else{die $$V{"all"} }		
	}else{
		($cChr2,$cPos2)=($$V{CHROM},$$V{"INFO"}{"END"});
	}
	$cLen=abs($cPos2-$cPos1);
	
	$$V{"INFO"}{"SVLEN"}=$cLen if !exists($$V{"INFO"}{"SVLEN"}) and $cChr1 eq $cChr2;
	$$V{"INFO"}{"SVLEN"}=abs($$V{"INFO"}{"SVLEN"}) if exists($$V{"INFO"}{"SVLEN"});
	
	$$V{INFO}{PAFG}=0;
	
	if(!exists($allowedChr{$cChr1}) ){print ${$V->print_line}; next}
	
	
	
	print ${$V->print_line} if $verbose;

	if($$V{"INFO"}{"TOOL"} eq "CNVnator"){ 
		print "--- $cChr1,$cPos1,$cPos2\n" if $verbose;
		$res = $tabix->query($cChr1,$cPos1,$cPos2); 	
	}else{ 
		$cPos1P=($cPos1+1000);
		$cPos1M=($cPos1-1000);
		$cPos1M=0 if $cPos1M<0;
		print "---MP $cChr1,$cPos1M,$cPos1P\n" if $verbose;
		$res = $tabix->query($cChr1,$cPos1M,$cPos1P); 
	}
	
	
	
	undef @allOvl;
	while(my $line = $tabix->read($res)){ #1       869477  870222  name=Lumpy_28;sample=FR05812673;SVTYPE=DEL;SVLEN=745;TOOL=Lumpy;SR=0;PE=8;DRF=0.12;IDD=.;       0       +       869477  870222  DEL     1       746     0
		
		($bChr,$bSt,$bEn,$bName)=split("\t",$line);

		undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=$cTMP2[1]; }
		
		($bChr1,$bPos1,$bChr2,$bPos2)=($bChr,($bSt+1),$h{"CHR2"},$bEn,);
		
		$bLen=($bEn-$bSt+1);
		
		print "$bChr1,$bPos1,$bChr2,$bPos2 \n"  if $verbose;
		#$tmp=<STDIN>; 
		
# 		next if $$V{"INFO"}{"SVTYPE"} ne $h{"SVTYPE"};

		
		##### check if breakpoints are within margins, else skip
		$brDist1=abs($cPos1-$bPos1);
		$brDist2=abs($cPos2-$bPos2);
		
# 		print "  test chromosome $bChr1 ne $cChr1 :: $line \n" if $verbose;
		next if $bChr1 ne $cChr1;

# 		print "  test chromosome $bChr2 ne $cChr2 :: $line \n" if $verbose;
		next if $bChr2 ne $cChr2;
		


		if (($h{"SVTYPE"} eq "BND" or $h{"SVTYPE"} eq "DEL" or $h{"SVTYPE"} eq "DUP" or $h{"SVTYPE"} eq "INV") 
			and  $h{"SVTYPE"} ne $$V{"INFO"}{"SVTYPE"}){
				next
		}

		
		# MCNV -> DEL/DUP
		next if $h{"SVTYPE"} eq "MCNV" and $$V{"INFO"}{"SVTYPE"} ne "DEL" and $$V{"INFO"}{"SVTYPE"} ne "DUP";
		
		if ($h{"SVTYPE"} eq "MCNV"){
			#print STDERR $h{AF}."\n";
			@AFs=split(",",$h{AF}); 
			$h{AF}=max(@AFs);
			#print STDERR $h{AF}."\n";
		}
		 
		
		# always pass CPX, CTX

		
		# for INS 
		next if $h{"SVTYPE"} eq "INS" and $$V{"INFO"}{"SVTYPE"} ne "BND";
		if ($h{"SVTYPE"} eq "INS"){
			$brDist2=abs($cPos2-$bPos1);
		}
		
		
		next if $brDist1>1000 or $brDist2>1000;
		
# 		if ($cChr1 eq $cChr2){ # calculate overlap if same chromosome
# 			$cOvl=overlap($cPos1,$cPos2,$bPos1,$bPos2);	
# 			
# 			$cOvlDiff1=abs($cLen-$cOvl);
# 			$cOvlDiff2=abs($bLen-$cOvl);
# 			$cOvlFrac1=sprintf("%.5f",($cOvl/$cLen));
# 			$cOvlFrac2=sprintf("%.5f",($cOvl/$bLen));	
# 		}
# 		
# 		# is CNV has to match
# 		if(exists($$V{"INFO"}{"CNV"}) and exists($h{"CNV"}) and $h{"CNV"} != $$V{"INFO"}{"CNV"}){next}
# 		
# 		if(exists($$V{"INFO"}{"CNV"}) and $$V{"INFO"}{"CNV"}==1 and $cChr1 eq $cChr2){ # if CNV
# 
# 			# at least 20% of the KDB variant has ot be covered
# 			if($cOvlFrac2<0.2){next} 
# 			# at least 70% of ClinSV variant has to be covered, or 90% if the length difference is >1000 bases
# 			if( ($cOvlDiff1>1000 and $cOvlFrac1<0.9) or $cOvlFrac1<0.7 ){ next }
# 			
# 		}else{ # balanced event
# 			next if $brDist1>1000 or $brDist2>1000;
# 						
# 		}

		#$cEvid=""; foreach $cType (("SR","PE","DRF","IDD","CNRD","CNP")){ $cEvid.="$cType=".$h{$cType}."," if exists($h{$cType}); }	chop($cEvid);
		
		push @allOvl, [$bChr1,$bPos1,$bChr1,$bPos2,($brDist1+$brDist2),$h{name}.",".$h{"SVTYPE"},$h{AC},$h{AF}];
		
		print "$cChr1:$cPos1\-$cPos2 KDB:$bChr1:$bPos1\-$bPos2 brDist1:$brDist1 brDist2:$brDist2, cOvl:$cOvl cOvlFrac1:$cOvlFrac1, cOvlFrac2:$cOvlFrac1 cOvlDiff1:$cOvlDiff1 cOvlDiff2:$cOvlDiff2|| $line\n" if $verbose;
		
		$$V{INFO}{PAFG}=roundA($h{AF});
		$presentInKDB++;
		break;
	}



	print ${$V->print_line};
	
	
# 	last if $c2++>100;
}

print STDERR "of all variants in VCF $allVariantsInVCF, $presentInKDB (".sprintf("%.1f",($presentInKDB/$allVariantsInVCF*100))."%) overlap with gnomAD\n";



sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}


sub roundA{ my ($cFreq)=@_;  if($cFreq!=1){ $RL=int(log10($cFreq)*(-1))+2; $cFreq=round($cFreq,$RL);  }  return $cFreq }
sub log10 { my $n = shift;  my $r=($n>0)? log($n)/log(10):0; return $r; }

sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}


















