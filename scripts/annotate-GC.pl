
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use Compress::Zlib;

$inVCF=$ARGV[0];
$Ref=$ARGV[1];

use My::Seq qw(GC_seq);
use My::VCF2;
 
############ input sorted vcf

my $VCFObj=My::VCF2->new($inVCF);
$VCFObj->add_header("##INFO=<ID=GC,Number=1,Type=Integer,Description=\"GC content of the variant\">"); 
$VCFObj->add_header("##INFO=<ID=CR,Number=1,Type=Float,Description=\" Size ratio of compressed (zip) vs. uncompressed reference sequence of the variant. Low complexity sequences have smaller compression ratios.\">"); 
print ${$VCFObj->print_header};

$lChr="";$cSeq="";
while (my $V=$VCFObj->next_line){	
	
	($cChr,$cSt,$cEn)=($$V{CHROM},$$V{POS},$$V{"INFO"}{"END"});		
	

	if (exists($$V{"INFO"}{"GC"})){ print ${$V->print_line}; $cC++; next}
	
	if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
		# 		next if exists($$V{"INFO"}{"SECONDARY"});
		if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			($cChr2,$cPos2)=($1,$2);
		}else{die "can not interprete ALT field".$$V{"all"} }
	}else{
		($cChr2)=($$V{CHROM});
	}
	if ($cChr ne $cChr2){ print ${$V->print_line}; $cC++; next}
	
	
	
	if($$V{"INFO"}{"CNV"} == 0){ print ${$V->print_line}; $cC++; next;}
	
	$cReg=($cSt<$cEn)? "$cChr:$cSt-$cEn":"$cChr:$cEn-$cSt";
# 	print STDERR "samtools faidx $Ref $cReg\n";
	$cSeq=""; open(IN,"samtools faidx $Ref $cReg |"); while(<IN>){next if $.==1; chomp; $cSeq.=$_}close(IN);
# 	print STDERR $cSeq."\n";
	
	
	$$V{"INFO"}{"CR"}=(length($cSeq)>0)? sprintf("%.3f",length(compress($cSeq))/length($cSeq)):(-1);
	
	$$V{"INFO"}{"GC"}=sprintf("%.0f",GC_seq(\$cSeq));
# 	print "".sprintf("%.0f",GC_seq(\$cSeq))." ".length($cSeq)."\n";
	print ${$V->print_line};
	
	$cC++;
	$lChr=$cChr;
	
}


# print STDERR "$cC\n";



























