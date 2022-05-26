
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";


use My::InterS qw(merge_segs);
use My::VCF2; 
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;

($inF,$inVCF)=@ARGV;

%allowedChr;

%DGV2LumpyType=("duplication","DUP","gain","DUP","deletion","DEL","loss","DEL","inversion","INV");

#   89612 deletion
#   69524 loss
#   26384 gain
#   25596 insertion
#   21360 duplication
#    8974 novel sequence insertion
#    5759 gain+loss
#    4156 mobile element insertion
#    2000 sequence alteration
#    1328 inversion
#     517 tandem duplication
#     466 complex


my $IntersObj=My::InterS->new;


open(IN1, "gzip -dc $inF | ") || die "$inF nicht gefunden";
while(<IN1>){ next if /^#/;
	chomp; @_=split("\t",$_);
	$allowedChr{$_[0]}++;
}close(IN1);


############################# annotate the VCF file




my $tabix = new Tabix(-data => "$inF",-index => "$inF.tbi");


#definition flanking seg-dups: max 20% of len at +- 20% of variant ends


my $VCFObj=My::VCF2->new($inVCF);

$VCFObj->add_header("##INFO=<ID=SEGD,Number=1,Type=String,Description=\"Overlapping segmental duplications published by Bailey JA et al. 2002. For best match: % variant coverage | % seg-dup coverage | identity | for all matching seg-dup's: count | merged % variant coverage\">");

print ${$VCFObj->print_header};

while (my $V=$VCFObj->next_line){
  	
  	($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{"INFO"}{"END"});
  	
	
	if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
		# 		next if exists($$V{"INFO"}{"SECONDARY"});
		if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			($cChr2,$cPos2)=($1,$2);
		}else{die "can not interprete ALT field".$$V{"all"} }
	}else{
		($cChr2)=($$V{CHROM});
	}
	
  	if( $$V{"INFO"}{"CNV"} == 0 or $cChr2 ne $cChr or !exists($allowedChr{$cChr}) ){ print ${$V->print_line}; next } # do not annotate copy neutral events
  	
	($cSt,$cEn)=($cSt<$cEn)? ($cSt,$cEn):($cEn,$cSt);
	$cLen=$cEn-$cSt+1;
# 	print "--- $cChr,$cSt,$cEn\n" ;
	my $res = $tabix->query($cChr,$cSt,$cEn);
	
	undef @allOvl;
	while(my $line = $tabix->read($res)){ #1        13144727        13150117        name=147_2;otherLoc=chr1:13047629-13052998;fracMatch=0.96386;   0       +       13144727        13150117        182,182,182     1       5390    0
		
		($bChr,$bSt,$bEn,$bName)=split("\t",$line);
		$bSt++;
		$bLen=($bEn-$bSt+1);
		
		undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=$cTMP2[1]; }
		$cOvl=overlap($cSt,$cEn,$bSt,$bEn);
		$cOvlFrac1=sprintf("%.5f",($cOvl/$cLen));
		$cOvlFrac2=sprintf("%.5f",($cOvl/$bLen));
		
		push @allOvl, [$bChr,$bSt,$bEn,$bLen,$cOvl,$cOvlFrac1,$cOvlFrac2,$h{fracMatch},($cOvlFrac1+$h{fracMatch}),$h{name},$h{otherLoc}];
		
# 		print "$cChr,$cSt,$cEn || $line\n";
	}
	@allOvl=sort {$$b[8] <=> $$a[8] } @allOvl;
	$ovlCount=scalar(@allOvl);
# 	foreach $cA (@allOvl){ 	print "Best: ".join("\t",@$cA)."\n"; }	
	if(scalar(@allOvl)>0){ $bestM= sprintf("%.0f",$allOvl[0][5]*100)."|".sprintf("%.0f",$allOvl[0][6]*100)."|".sprintf("%.1f",$allOvl[0][7]*100)."|".scalar(@allOvl);}else{$bestM="0|0|0|0"}
	
	merge_segs(\@allOvl);
	$ovlSum=0; foreach $cA (@allOvl){ $ovlSum+=overlap($cSt,$cEn,$$cA[1],$$cA[2]); }
	$ovlPc=sprintf("%.0f",($ovlSum/$cLen)*100);
	
	################ check flanking
	
	#$hasFlankingDelDups=checkFlaning($cChr,$cSt,$cEn,$cLen);
	
	if($ovlCount>0){ $$V{INFO}{SEGD}=$bestM."|$ovlPc";}
	
	print ${$V->print_line};
	
	
# 	last if $c2++>100;
}





sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}















































