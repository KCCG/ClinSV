
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
use Getopt::Long;
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;

$concDist=1000;
$cnvPcOvlOfDb=20; 
$cnvPcOvlOfVp=90;
$cnvPcOvlOfVp2=70;

GetOptions ("dgv=s"  => \$inDGV,
			"toAnnot=s"  => \$toAnnot,
			"distance=f" => \$concDist,
			"ref=s" => \$inRef
)  or die("Error in command line arguments\n");

$verbose=0;
print STDERR "# Summed distance to break-point pair at which two SV are still concordant: $concDist\n";
print STDERR "# toAnnot $toAnnot\n";


#############################
############################# 	1. prepare 1k genome DGV flie
#############################


# 			466 complex
# 			517 tandem duplication
# 		   1428 inversion
# 		   2000 sequence alteration
# 		   4156 mobile element insertion
# 		   6902 gain+loss
# 		   8974 novel sequence insertion
# 		  21360 duplication
# 		  25610 insertion
# 		  46452 gain
# 		 101257 deletion
# 		 121527 loss

%type2col=("DEL","#c34343","loss","#c34343","deletion","#c34343",
"DUP","#0375b6","gain","#0375b6","duplication","#0375b6","tandem duplication","#0375b6",
"INV","#a968b5","inversion","#a968b5",
"TRA","#f8a748","insertion","#f8a748","novel sequence insertion","#f8a748","mobile element insertion","#f8a748",
"complex","#8d4507","sequence alteration","#8d4507","gain+loss","#8d4507");




# read in the allowed Chr for tabix, else error
%allowedChr;
print STDERR "# read in the allowed Chr for DGV tabix file, else error if try to access chr that does not exist... ";
open(IN1, "gzip -dc $inDGV | ") || die " $inDGV not found";
while(<IN1>){ chomp; @_=split("\t",$_);
	$allowedChr{$_[0]}++;
}
close(IN1);
print STDERR "  ".join(",",sort keys %allowedChr)." allowed chrs. \n";

$tabix = new Tabix(-data => "$inDGV",-index => "$inDGV.tbi");


# read in chr 2 len
open(IN1, "<$inRef.fai") || die " $inDGV not found";
while(<IN1>){ chomp; @_=split(/[\t| ]/,$_);
	$chr2len{$_[0]}=$_[1];
}
close(IN1);






#############################
############################# 2. loop through B
#############################

print STDERR "# start annotating... \n";
%chrNotExist;
$concV=0;
our %Aconc=();
if ($toAnnot =~ /vcf$|vcf.gz$/){	
      
    $V=My::VCF2->new($toAnnot);
    
    $outVCF=$toAnnot; $outVCF=~ s/vcf[.]gz$|vcf$//g; $outVCF.="DGV.vcf";
    open(VCF, ">$outVCF") || die "can not open $outVCF";
    print STDERR "# out vcf $outVCF \n";

	$V->add_header("##INFO=<ID=PAFDGV,Number=.,Type=Float,Description=\" DGV population variant allele frequency of best hit \">");
	$V->add_header("##INFO=<ID=DGVD,Number=.,Type=String,Description=\" Top 3 DGV hits VAF, study sample size, sum reciprocal % overlap, DGV ID, publication, detection method \">");
	print VCF ${$V->print_header};

	while (my $V=$V->next_line){
		
		$allVariantsinA++;
	    my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV)=obtainCoordsFromVCF($V);
	    push @{$stat{countB}}, $cLen;
# 		print STDERR ${$V->print_line} if $verbose;
		
		
		
		
		if(!exists($allowedChr{$cChr1}) ){print VCF ${$V->print_line}; next}
		
		#print "$cChr1 ne $cChr2 \n";
		
		$$V{INFO}{PAFDGV}=(-1);
		$resConc=conCheck($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV);
		
		print VCF ${$V->print_line};

	}
}



#############################
############################# 3. print summary
#############################

print STDERR "$concV of $allVariantsinA variants could be annotated with DGV\n";
print STDERR "# ".join(",",sort keys %chrNotExist)." skipped chrs. \n";

#############################
############################# Functions
#############################


sub	conCheck {
	
	my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV)=@_;
	undef %hits;
	
	return if $cType eq "BND";
	
	$concFound=0;		
	
	($xChr1,$xPos1,$xChr2,$xPos2,$xLen,$xType,$xCNV)=($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV);
	
    # 2->5, 1 -> 10, 3->3
    
    if (!exists($allowedChr{$xChr1})){  $chrNotExist{$xChr1}++; return;}


	$cMin=($xPos1-($xLen*100/$cnvPcOvlOfDb));
	$cMax=($xPos2+($xLen*100/$cnvPcOvlOfDb));
	$cMin=1 if $cMin<0;
	$cMax=$chr2len{$xChr1} if $cMax>$chr2len{$xChr1};
	
	print "    ---vcf $xChr1:$xPos1\-$xPos2 $xType $cMin,$cMax\n" if $verbose;
	
	$res = $tabix->query($xChr1,$cMin,$cMax);
	
	$passCat="tabix2";
	while(my $line = $tabix->read($res)){ #1       869477  870222  name=Lumpy_28;sample=FR05812673;SVTYPE=DEL;SVLEN=745;TOOL=Lumpy;SR=0;PE=8;DRF=0.12;IDD=.;       0       +       869477  870222  DEL     1       746     0
	
		#1	1443537	1445747	type=copy_number_loss;typeid=SO:0001743;varid=esv3818131;count=108;freq=0.043
		chomp($line); ($bChr1,$bPos1,$bPos2,$bType)=split("\t",$line); $bPos1++; $bLen=($bPos2-$bPos1+1);
		undef %h; foreach $cT (split(";",$bType)){ @t2=split("=",$cT); $h{$t2[0]}=$t2[1];  }
		
# 			466 complex
# 			517 tandem duplication
# 		   1428 inversion
# 		   2000 sequence alteration
# 		   4156 mobile element insertion
# 		   6902 gain+loss
# 		   8974 novel sequence insertion
# 		  21360 duplication
# 		  25610 insertion
# 		  46452 gain
# 		 101257 deletion
# 		 121527 loss

		if($h{type} eq "deletion" or $h{type} eq "loss"){ $bType="DEL";  }
		elsif($h{type} eq "gain" or $h{type} eq "duplication" ){ $bType="DUP";}
		elsif($h{type} eq "mobile_element_insertion" or $h{type} eq "insertion"){ $bType="BND";}
		elsif($h{type} eq "inversion" ){ $bType="INV";}
		else{next}
		
# 		print "      tabix2 --- $bChr1:$bPos1\-$bPos2,$bLen,$bType\n";
# 		
		$passCat="CNV"; next if ($h{type} eq "loss" or $h{type} eq "gain") and $xCNV==0; # if loss or gain, the vcf variant can not be copy neutral

		$passCat="type"; next if $xType ne $bType; # type DEL, INV, .. has to match
		
		$passCat="chr"; next if $bChr1 ne $xChr1; # chr have to match
		
		$cDist1=abs($bPos1-$xPos1); $cDist2=abs($bPos2-$xPos2); $cSumDist=($cDist1+$cDist2);
		
		$cOvl=overlap($bPos1,$bPos2,$xPos1,$xPos2);
		$xOvlDiff=abs($xLen-$cOvl);
		$bOvlDiff=abs($bLen-$cOvl);
		$xOvlPc=sprintf("%.1f",($cOvl/$xLen*100));
		$bOvlPc=sprintf("%.1f",($cOvl/$bLen*100));
		$cSumPc=($xOvlPc+$bOvlPc);
		
# 			print "$xOvlPc $bOvlPc \n";
		$passCat="ovl";	
		if($cType eq "INV" or $xCNV==0){ next if $cDist1>1000 or $cDist2>1000; } # if ClinSV balanced
		else{ # CNVs
			# at least 20% of the DGV variant has ot be covered
			if($bOvlPc<$cnvPcOvlOfDb){next} 
			# at least 70% of ClinSV variant has to be covered, or 90% if the length difference is >1000 bases
			if( ($xOvlDiff>1000 and $xOvlPc<$cnvPcOvlOfVp) or $xOvlPc<$cnvPcOvlOfVp2 ){ next }
		}
		 # pos not further away than concDist
# 			print "$xOvlPc $bOvlPc  -- pass\n";
					
		$concFound++;
		
		#id=variantaccession;type=$varianttype;reference=$reference;pubmedid=$pubmedid;method=$method;freq=$cFreq;samplesize=$samplesize
		$hits{$h{id}}={%h};
		$hits{$h{id}}{cSumPc}=$cSumPc;
# 		print "concordance\toverlap\t$cSumPc\t$bChr1:$bPos1-$bPos2,$bType\t$xChr1:$xPos1-$xPos2,$xType\t$bOvlPc,$xOvlPc\t$bLen,$cLen\t$line\n" if $verbose;
	}
	
	# assign freq of best hit
	undef @oDGVD;
	foreach $cID (sort { $hits{$b}{cSumPc} <=> $hits{$a}{cSumPc} } keys %hits){
		$$V{INFO}{PAFDGV}=$hits{$cID}{freq} if scalar(@oDGVD)==0;
		push @oDGVD, join(",",($hits{$cID}{freq},$hits{$cID}{samplesize},round($hits{$cID}{cSumPc},0)."pc",$hits{$cID}{id},$hits{$cID}{reference},$hits{$cID}{"method"}));
		last if scalar(@oDGVD)>=3;
	}
	$$V{INFO}{DGVD}=join("|",@oDGVD) if @oDGVD>0;
	
	$concV++ if $concFound>0;
# 	print " passCat: $passCat \n" if $verbose;

}


sub obtainCoordsFromVCF{
	
	($X)=@_;
	($aChr1,$aPos1)=($$X{CHROM},$$X{POS});
	
	if ($$X{"ALT"} =~ /del|DEL/){$aType="DEL"}
	elsif ($$X{"ALT"} =~ /dup|DUP/){$aType="DUP"}
	elsif ($$X{"ALT"} =~ /inv|INV/){$aType="INV"}
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
	$aCNV=$$X{"INFO"}{"CNV"};
	return((($aChr1,$aPos1,$aChr2,$aPos2,$aLen,$aType,$aCNV)));
}


sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}



sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}

sub log10 { my $n = shift;  return log($n)/log(10); }























