
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
%fh;
$outStem=shift(@ARGV);

print IGV "#gffTags\n";
%PE_type2col=("DEL","#c34343","DUP","#0375b6","INV","#a968b5","TRA","#f8a748","BND","#f8a748");

%exclInfoTag=("EVTYPE",1,"IMPRECISE",1,"SEGD",1,"GDESC",1,"sample",1,"KDBB",1,"KDBD",1,"GT",1);

foreach $inVCF (@ARGV){
my $VCFObj=My::VCF2->new($inVCF); 
	@samples=@{$$VCFObj{samples}};
	
	while (my $V=$VCFObj->next_line){
		
		
		if ($$V{INFO}{SVTYPE} eq "BND"){
# 			next if exists($$V{INFO}{"SECONDARY"});
			if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
				
				$$V{INFO}{CHROM2}=$1;
				$$V{INFO}{POS2}=$2;
				if($$V{CHROM} eq $1){
					$$V{INFO}{SVLEN}=abs($$V{POS}-$$V{INFO}{POS2});
				}else{
					$$V{INFO}{SVLEN}=0;
				}
			}else{die $$V{"all"} }
		}else{
			$$V{INFO}{CHROM2}=$$V{CHROM};
			$$V{INFO}{POS2}=$$V{INFO}{END};
			
			if(exists($$V{INFO}{SVLEN})){
				$$V{INFO}{SVLEN}=abs($$V{INFO}{SVLEN});
			}else{
				$$V{INFO}{SVLEN}=abs($$V{POS}-$$V{INFO}{END});
			}	
		}

		
		$cCol="#5d5553"; $cCol=$PE_type2col{$$V{INFO}{SVTYPE}} if exists ($PE_type2col{$$V{INFO}{SVTYPE}});

		foreach $cSample ( @{$$V{samples}} ){
			next if $$V{GTS}{GT}{$cSample} eq "0/0" or $$V{GTS}{GT}{$cSample} eq "./."; 
			next if $$V{GTS}{FT}{$cSample} eq "FAIL";
			
			if(!exists($fh{KDB}{$cSample})){ open($fh{KDB}{$cSample}," > $outStem.$cSample.bed") || die "! $outStem.$cSample.bed" ; print {$fh{KDB}{$cSample}} "#gffTags\n"; }
			if(!exists($fh{IGV}{$cSample})){ open($fh{IGV}{$cSample}," > $outStem.$cSample.igv.bed") || die "! $outStem.$cSample.igv.bed" ; print {$fh{IGV}{$cSample}} "#gffTags\n"; }
			
			$gffTags="name=".$$V{ID}.",$cSample;";
			
			if($$V{INFO}{SVTYPE} eq "BND"){
				$gffTags.="SV=".$$V{"INFO"}{SVTYPE}; 
				$gffTags.=",LEN:".$$V{"INFO"}{SVLEN} if exists($$V{"INFO"}{SVLEN}); 
				$gffTags.=",LOC:".$$V{"INFO"}{LOCATION} if exists($$V{"INFO"}{LOCATION});
				$gffTags.=";";
			}else{
				$gffTags.="SV=".$$V{"INFO"}{SVTYPE}.",LEN:".$$V{"INFO"}{SVLEN};
				$gffTags.=",LOC:".$$V{"INFO"}{LOCATION} if exists($$V{"INFO"}{LOCATION});
				$gffTags.=";";
			}
			
			$gffTags.="TOOL=".$$V{"INFO"}{TOOL}.",FT:$$V{GTS}{FT}{$cSample},RARE:$$V{GTS}{RARE}{$cSample};";
			
			undef %hGTS;undef @aGTS;
			push @aGTS, ("PE","SR","MQBP");
			if($$V{"INFO"}{TOOL} =~ /Lumpy/){ push @aGTS, ("DBP","IDD") }
			if($$V{"INFO"}{TOOL} =~ /CNVnator/i){ push @aGTS, ("CNP") }
			
			$gffTags.="EVID="; 
			foreach $cT (@aGTS){ 
				next if $hGTS{$cT}++>0;
				$gffTags.="$cT:".$$V{GTS}{$cT}{$cSample}."," if exists($$V{GTS}{$cT}) and exists($$V{GTS}{$cT}{$cSample});  
			} chop($gffTags);
			$gffTags.=",CNV:".$$V{"INFO"}{CNV} if exists($$V{"INFO"}{CNV}); 
			$gffTags.=",GC:".$$V{"INFO"}{GC} if exists($$V{"INFO"}{GC}); 
			$gffTags.=",MQ:".$$V{"INFO"}{MQ} if exists($$V{"INFO"}{MQ}); 
			$gffTags.=";";
			
			@aGTS=("DRF","DRA");
			$gffTags.="DEPTH="; 
			foreach $cT (@aGTS){ 
				next if $hGTS{$cT}++>0;
				$gffTags.="$cT:".$$V{GTS}{$cT}{$cSample}."," if exists($$V{GTS}{$cT}) and exists($$V{GTS}{$cT}{$cSample});  
			} chop($gffTags);
			$gffTags.=";";
							
			$gffTags.="PAFV=".$$V{"INFO"}{PAFV}.";" if exists($$V{"INFO"}{PAFV});  
			$gffTags.="PAFSU=".$$V{GTS}{PAFSU}{$cSample}.";" if exists($$V{GTS}{PAFSU}{$cSample});  
			$gffTags.="PAFDRA=".$$V{GTS}{PAFDRA}{$cSample}.";" if exists($$V{GTS}{PAFDRA}{$cSample});  
			$gffTags.="PAF1KG=".$$V{"INFO"}{PAF1KG}.";" if exists($$V{"INFO"}{PAF1KG});  
			$gffTags.="PAFG=".$$V{"INFO"}{PAFG}.";" if exists($$V{"INFO"}{PAFG});  
			$gffTags.="DGVD=".$$V{"INFO"}{DGVD}.";" if exists($$V{"INFO"}{DGVD});
			
			$gffTags.="ENS=".$$V{"INFO"}{NUMG}; 
			if (exists($$V{"INFO"}{GENES})){ $GENES=$$V{"INFO"}{GENES}; $GENES=~ s/[\" ]//g; $gffTags.=",".$GENES  }
			$gffTags.=";";
				
				
			########## print IGV, for IGV also pring BND as one segment if close	
			if ($$V{INFO}{CHROM2} ne $$V{CHROM} or ($$V{INFO}{SVLEN}>1000000 and $$V{"INFO"}{CNV}==0) or 
			($$V{INFO}{SVLEN}>1000000 and $$V{"INFO"}{CNV}==1 and $$V{"INFO"}{SVTYPE} eq "BND")
			){ # split breakpoint over two lines for igv
				
				($cSt,$cEn)=(($$V{POS}-1),$$V{POS});
				print {$fh{IGV}{$cSample}}  join("\t",($$V{CHROM},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
				
				($cSt,$cEn)=(($$V{INFO}{POS2}-1),$$V{INFO}{POS2});
				print {$fh{IGV}{$cSample}}  join("\t",($$V{INFO}{CHROM2},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
				
			}else{
				($cSt,$cEn)=(($$V{POS}-1),$$V{INFO}{POS2});
				print  {$fh{IGV}{$cSample}}  join("\t",($$V{CHROM},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
			}
			
			
			$gffTags="name=".$$V{ID}.";sample=$cSample;";
			$gffTags.="ALT=".$$V{"ALT"}.";" if $$V{INFO}{SVTYPE} eq "BND";
			foreach $cT (sort keys %{$$V{INFO}}){ $gffTags.="$cT=".$$V{INFO}{$cT}.";" if exists($$V{INFO}{$cT}); }
			foreach $cT (sort keys %{$$V{GTS}}){  $gffTags.="$cT=".$$V{GTS}{$cT}{$cSample}.";" if exists($$V{GTS}{$cT}) and exists($$V{GTS}{$cT}{$cSample}); }
			
			######## print KDB
			if ($$V{INFO}{SVTYPE} eq "BND"){
			
				($cSt,$cEn)=(($$V{POS}-1),$$V{POS});
				print  {$fh{KDB}{$cSample}}  join("\t",($$V{CHROM},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
				
				($cSt,$cEn)=(($$V{INFO}{POS2}-1),$$V{INFO}{POS2});
				print  {$fh{KDB}{$cSample}} join("\t",($$V{INFO}{CHROM2},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
				
			}else{
			
				($cSt,$cEn)=(($$V{POS}-1),$$V{INFO}{POS2});
				print  {$fh{KDB}{$cSample}}  join("\t",($$V{CHROM},$cSt,$cEn,$gffTags,0,"+",$cSt,$cEn,$cCol,1,$cLen,"0"))."\n";
			
			}
			
			
		}
		
	}

 
 
}

foreach $cSample ( @samples ){
close($fh{IGV}{$cSample});
print STDERR `sort_bgzip $outStem.$cSample.igv.bed`;
print STDERR `rm $outStem.$cSample.igv.bed`;
print STDERR `tabix -f -p bed $outStem.$cSample.igv.bed.gz`;
close($fh{KDB}{$cSample});
print STDERR `sort_bgzip $outStem.$cSample.bed`;
print STDERR `rm $outStem.$cSample.bed`;
print STDERR `tabix -f -p bed $outStem.$cSample.bed.gz`;
}



















































