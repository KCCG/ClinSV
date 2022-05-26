
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use Getopt::Long;

GetOptions ("force" => \$force,
			"infoFields=s" => \$infoFields,
			"gtsFields=s" => \$gtsFields,
		   "PASS"  => \$PASS,
		"outStem=s"  => \$outStem,
		"input=s"  => \$inVCF
)  or die("Error in command line arguments\n");

# print STDERR "inVCF: $inVCF\n";
@AinfoFields=split(",",$infoFields);
@AgtsFields=split(",",$gtsFields);

use My::VCF2; 
my $V=My::VCF2->new($inVCF); 

$aSamples=\@{$$V{samples}};

print STDERR "samples: ".join(",",@$aSamples)."\n";
print STDERR "outStem: $outStem\n";

%PE_type2col=("DEL","#c34343","DUP","#0375b6","INV","#a968b5","TRA","#f8a748");

%fh;

while (my $V=$V->next_line){
 ($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{"INFO"}{"END"}); $cLen=($cEn-$cSt+1);
 next if $$V{"INFO"}{"SVTYPE"} eq "BND";
 next if abs($$V{"INFO"}{"SVLEN"})>1000000;# and $$V{"INFO"}{"TOOL"} eq "Lumpy";
 $cCol="#5d5553"; $cCol=$PE_type2col{$$V{"INFO"}{"SVTYPE"}} if exists ($PE_type2col{$$V{"INFO"}{"SVTYPE"}});
 
 
 foreach $cSample ( @$aSamples ){
  
#   next if $$V{GTS}{GT}{$cSample} eq "0/0" or $$V{GTS}{GT}{$cSample} eq "./."; 
  next if $$V{GTS}{FT}{$cSample} ne "PASS";
#   print STDERR "$cChr,$cSt,$cEn".$$V{"INFO"}{"TOOL"}."\n";
  $gffTags="name=".$$V{ID}.";sample=$cSample;";
  foreach $cT (@AinfoFields){ $gffTags.="$cT=".$$V{"INFO"}{$cT}.";" if exists($$V{"INFO"}{$cT}); }
  foreach $cT (@AgtsFields){ $gffTags.="$cT=".$$V{GTS}{$cT}{$cSample}.";" if exists($$V{GTS}{$cT}{$cSample}); }
#   print STDERR "$outStem.$cSample.bed.gz\n";
  if(!exists($fh{$cSample})){ open($fh{$cSample},"> $outStem.$cSample.bed") || die "! $outStem.$cSample.bed.gz" ; print {$fh{$cSample}} "#gffTags\n"; }
  print {$fh{$cSample}} join("\t",($cChr,($cSt-1),$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
 }
}

foreach $cSample ( @$aSamples ){
close($fh{$cSample});
print STDERR `sort_bgzip $outStem.$cSample.bed`;
print STDERR `tabix -f -p bed $outStem.$cSample.bed.gz`;
}












