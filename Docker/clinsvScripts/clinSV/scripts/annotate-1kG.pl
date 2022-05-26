
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
use Getopt::Long;
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use File::Basename;

$concDist=1000;

$cnvPcOvlOfDb=20;
$cnvPcOvlOfVp=90;
$cnvPcOvlOfVp2=70;

GetOptions ("dgva=s"  => \$inDGVa,
			"toAnnot=s"  => \$toAnnot,
			"distance=f" => \$concDist
)  or die("Error in command line arguments\n");

$verbose=1;
print STDERR "# Summed distance to break-point pair at which two SV are still concordant: $concDist\n";
print STDERR "# toAnnot $toAnnot\n";


#############################
############################# 	1. prepare 1k genome DGVa flie
#############################


%type2col=("DEL","#c34343","copy_number_loss","#c34343","deletion","#c34343",
"DUP","#0375b6","copy_number_gain","#0375b6","duplication","#0375b6","tandem duplication","#0375b6",
"INV","#a968b5","inversion","#a968b5",
"TRA","#f8a748","insertion","#f8a748","mobile_element_insertion","#f8a748","insertion","#f8a748",
"complex","#8d4507","sequence alteration","#8d4507","indel","#8d4507");


$outDGVa=$inDGVa;# $outDGVa=~ s/[.]gvf[.]gz/.bed/g;


# read in the allowed Chr for tabix, else error
%allowedChr;
print STDERR "# read in the allowed Chr for DGVa tabix file, else error if try to access chr that does not exist... ";
open(IN1, "gzip -dc $outDGVa | ") || die " $outDGVa not found";
while(<IN1>){ chomp; @_=split("\t",$_);
	$allowedChr{$_[0]}++;
}
close(IN1);
print STDERR "  ".scalar(keys %allowedChr)." allowed chrs. \n";

$tabix = new Tabix(-data => "$outDGVa",-index => "$outDGVa.tbi");


#############################
############################# 2. loop through B
#############################

print STDERR "# start annotating... \n";

$concV=0;
our %Aconc=();
if ($toAnnot =~ /vcf$|vcf.gz$/){	
      
    $V=My::VCF2->new($toAnnot);
    
    $outVCF=$toAnnot; $outVCF=~ s/vcf[.]gz$|vcf$//g; $outVCF.="1kG.vcf";
    open(VCF, ">$outVCF") || die "can not open $outVCF";
    print STDERR "# out vcf $outVCF \n";

	$V->add_header("##INFO=<ID=PAF1KG,Number=.,Type=Float,Description=\" Variant frequency of 1k genome project consisting of 2504 genomes \">");
	print VCF ${$V->print_header};

	while (my $V=$V->next_line){
		
		$allVariantsinA++;
	    my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV)=obtainCoordsFromVCF($V);
	    push @{$stat{countB}}, $cLen;
# 		print STDERR ${$V->print_line} if $verbose;
		
		if(!exists($allowedChr{$cChr1}) ){print VCF ${$V->print_line}; next}
		

		$$V{INFO}{PAF1KG}=conCheck($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV);
		
		print VCF ${$V->print_line};

	}
}elsif($toAnnot =~ /bed$/){	

	$outVCF=$toAnnot; $outVCF=~ s/bed[.]gz$|bed$//g; $outVCF.="1kG.bed";
    open(OUTB, ">$outVCF") || die "can not open $outVCF";
	open(IN1, "<$toAnnot") || die " nicht gefunden";
	
	while(<IN1>){ chomp; 
	
		($cChr1,$cPos1,$cPos2,$cType)=split("\t",$_);
		$cPos1++;
		$cLen=$cPos2-$cPos1+1;
	
		$resConc=conCheck($cChr1,$cPos1,$cChr1,$cPos2,$cLen,$cType,1);	
		$cPos1--;
		print OUTB join("\t",($cChr1,$cPos1,$cPos2,$cType,$resConc))."\n";		
	}
	
	close(IN1);
	close(OUTB);
}

#############################
############################# 3. print summary
#############################

print STDERR "$concV of $allVariantsinA could be annotated with 1kG\n";



#############################
############################# Functions
#############################

sub	conCheck {
	
	my ($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV)=@_;
	undef %hits;
	
	$concFound=0;		
	return (-1) if $cType eq "BND";
	
	($xChr1,$xPos1,$xChr2,$xPos2,$xLen,$xType,$xCNV)=($cChr1,$cPos1,$cChr2,$cPos2,$cLen,$cType,$cCNV);
# 	print "    ---vcf $xChr1:$xPos1\-$xPos2 $xType $xLen\n";# if $verbose;
    
    # 2->5, 1 -> 10, 3->3
    
    $xPos1Min=($xPos1-($xLen*100/$cnvPcOvlOfDb));
    $xPos2Max=($xPos2+($xLen*100/$cnvPcOvlOfDb));
    $xPos1Min=1 if $xPos1Min<0;
    
	$res = $tabix->query($xChr1,$xPos1Min,$xPos2Max);
	
	$passCat="tabix2";
	if(defined $res->get){
	while(my $line = $tabix->read($res)){ #1       869477  870222  name=Lumpy_28;sample=FR05812673;SVTYPE=DEL;SVLEN=745;TOOL=Lumpy;SR=0;PE=8;DRF=0.12;IDD=.;       0       +       869477  870222  DEL     1       746     0
	
		#1	1443537	1445747	type=copy_number_loss;typeid=SO:0001743;varid=esv3818131;count=108;freq=0.043
		chomp($line); ($bChr1,$bPos1,$bPos2,$bType)=split("\t",$line); $bPos1++; $bLen=($bPos2-$bPos1+1);
		undef %h; foreach $cT (split(";",$bType)){ @t2=split("=",$cT); $h{$t2[0]}=$t2[1];  }
		
# 			168 insertion
# 			786 inversion
# 		   9132 indel
# 		   9414 copy_number_gain
# 		  16631 mobile_element_insertion
# 		  35939 copy_number_loss

		if($h{type} eq "copy_number_loss" or $h{type} eq "indel"){ $bType="DEL";}
		elsif($h{type} eq "copy_number_gain" ){ $bType="DUP";}
		elsif($h{type} eq "mobile_element_insertion" or $h{type} eq "insertion"){ $bType="BND";}
		elsif($h{type} eq "inversion" or $h{type} eq "insertion"){ $bType="INV";}
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
		$bOvlPc=sprintf("%.1f",($cOvl/$bLen*100)); # % of DGV variant
		$cSumPc=($xOvlPc+$bOvlPc);
		
# 			print "$xOvlPc $bOvlPc, $xOvlDiff,$bOvlDiff\n";
		$passCat="ovl";	
		if($cType eq "INV" or $xCNV==0){ next if $cDist1>1000 or $cDist2>1000; } # if ClinSV balanced
		else{ # CNVs
			# at least 20% of the DGV variant has ot be covered
			if($bOvlPc<$cnvPcOvlOfDb){next} 
			# at least 70% of ClinSV variant has to be covered, or 90% if the length difference is >1000 bases
			if( ($xOvlDiff>1000 and $xOvlPc<$cnvPcOvlOfVp) or $xOvlPc<$cnvPcOvlOfVp2 ){ next }
		}
		
		 # pos not further away than concDist
# 			print "$xOvlPc $bOvlPc  $cnvPcOvlOfDb $cnvPcOvlOfVp -- pass\n";
					
		$concFound++;
		$hits{$cSumPc}=$h{freq};
# 		print "concordance\toverlap\t$cSumPc\t$bChr1:$bPos1-$bPos2,$bType\t$xChr1:$xPos1-$xPos2,$xType\t$bOvlPc,$xOvlPc\t$bLen,$cLen\t$line\n" if $verbose;
	}
	}
	# assign the smallest freq
	foreach $cSumPc (sort { $hits{$a} <=> $hits{$b} } keys %hits){
		$concV++ if $concFound>0;
		return $hits{$cSumPc};
		last;
	}
	
	return (0);
	
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
	$aCNV=$$X{"INFO"}{CNV};
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

































