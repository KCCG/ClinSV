
# Author: Andre Minoche, a.minoche@garvan.org.au



use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";
use Excel::Writer::XLSX;
use My::VCF2; 
use File::Basename;

$inputStem=shift(@ARGV); 
$projectDir=shift(@ARGV); 
$S_Phen=shift(@ARGV);
$IGVSessionDir=shift(@ARGV); 
$inRef=shift(@ARGV); 

if (basename($inRef) =~ /38/){
	$chrPf='chr';
	$Y_chr="chrY";
	$X_chr="chrX";
	$M_chr="chrM";
}else{
	$chrPf='';
	$Y_chr="Y";
	$X_chr="X";
	$M_chr="MT";
}
###### read the ensemble to PHEN database 
open(IN1, "<$S_Phen") || die " $S_Phen not found";
while(<IN1>){ chomp; @_=split("\t",$_);
	next if $.==1;
	$ph_gene2type2desc{$_[4]}{$_[6]}{$_[5]}++;
	$ph_gene2mim{$_[4]}{$_[1]}++ if $_[1] ne "";	
	
}close(IN1);


foreach $cGene ( keys %ph_gene2type2desc ){
	if (  exists($ph_gene2mim{$cGene})  ){
		$ph_gene2phen{$cGene}="MIM - ".join("; ",keys %{$ph_gene2mim{$cGene}});
	}elsif (  exists($ph_gene2type2desc{$cGene}{"MIM disease"})  ){
		$ph_gene2phen{$cGene}="MIM - ".join("; ",keys %{$ph_gene2type2desc{$cGene}{"MIM disease"}});
	}elsif (  exists($ph_gene2type2desc{$cGene}{"DDG2P"})  ){
		$ph_gene2phen{$cGene}="DDG2P - ".join("; ",keys %{$ph_gene2type2desc{$cGene}{"DDG2P"}});
	}elsif (  exists($ph_gene2type2desc{$cGene}{"Orphanet"})  ){
		$ph_gene2phen{$cGene}="Orphanet - ".join("; ",keys %{$ph_gene2type2desc{$cGene}{"Orphanet"}});
	}
	
# 	print STDERR "$cGene ".$ph_gene2phen{$cGene}."\n";
}

undef %ph_gene2type2desc; undef %ph_gene2mim;



######## if ped fiel is present extract Affected onAffected
%sa2info;
if(-e "$projectDir/sampleInfo.ped"){
	$ped=1;
	print STDERR "reading from ped $projectDir/sampleInfo.ped ...\n";
	open(PED, "<$projectDir/sampleInfo.ped") || die "can not write to $projectDir/sampleInfo.ped\n";
		while(<PED>){
			next if /^#/;
			chomp; @t=split("\t",$_);
			$sa2Aff{$t[1]}++ if $t[5] == 2;
			$sa2nonAff{$t[1]}++ if $t[5] == 1;
			$sa2info{$t[1]}{"affected"}=$t[5];	
			$sa2info{$t[1]}{"family"}=$t[0];	
			$infoFields=0; if(@t>5){  map { $infoFields++;    $sa2info{$t[1]}{"pedInfo$infoFields"}=$t[$_];   } 6..$#t;  } # print STDERR "$t[1] pedInfo$infoFields $t[$_] \n";
# 			
		}
	close(PED);
}


if(-e "$projectDir/testGene.ids" or -e "$projectDir/candGene.ids"){
	$geneL=1;
	
	if(-e "$projectDir/candGene.ids"){ 
		print STDERR "reading gene list $projectDir/candGene.ids ...\n";
		open(GLIST, "<$projectDir/candGene.ids") || die "can not read to $projectDir/candGene.ids\n"; }
	elsif(-e "$projectDir/testGene.ids"){ 
		print STDERR "reading gene list $projectDir/testGene.ids ...\n";
		open(GLIST, "<$projectDir/testGene.ids") || die "can not write to $projectDir/testGene.ids\n"; }
	
		while(<GLIST>){		
			chomp; @t=split("\t",$_);
			$geneList{$t[0]}++;
# 			print "$t[0]\n";
		}
	close(GLIST);
}


print scalar(keys %geneList)." genes found\n";



print STDERR "read input vcf and prioritize variants: $inputStem.vcf ...\n";
	
open(VCF, ">$inputStem.prioritized.vcf") || die "can not write to $inputStem.vcf\n";
open(TXT, ">$inputStem.prioritized.txt") || die "can not write to $inputStem.txt\n";
open(RAREP, ">$inputStem.prioritized.RARE_PASS_GENE.txt") || die "can not write to $inputStem.txt\n";
open(RAREPV, ">$inputStem.prioritized.RARE_PASS_GENE.vcf") || die "can not write to $inputStem.txt\n";
open(PASS, ">$inputStem.prioritized.PASS.vcf") || die "can not write to $inputStem.txt\n";

$xl{wb}{RGfull} = Excel::Writer::XLSX->new("$inputStem.prioritized.RARE_PASS_GENE.xlsx");
$xl{wb}{RGlight} = Excel::Writer::XLSX->new("$inputStem.prioritized.RARE_PASS_GENE.light.xlsx");
$xl{wb}{full} = Excel::Writer::XLSX->new("$inputStem.prioritized.xlsx");

map { $xl{ws}{$_} = $xl{wb}{$_}->add_worksheet(); } (("RGfull","RGlight","full"));

my $V=My::VCF2->new($inputStem.".vcf"); 
$V->add_header("##FORMAT=<ID=RARE,Number=1,Type=Integer,Description=\"Is variant rare? 1=yes, 0=no Rare means that PAF, PAFSUF, PAFDRA and PAF1KG are smaller or equal 1%\">");

$V->add_header("##FILTER=<ID=LOW,Description=\"Variant does not pass in any sample, see FT\">");
$V->add_header("##FILTER=<ID=PASS,Description=\"Variant passing in at least one of the samples, see FT\">");


@columns=("SAMPLE","ID","FT","RARE","SU","PAFSU","PE","SR","DRF","DRA","PAFDRA","PCSD","GT","MQBP","CNV","IGV","GOTO","LOCATION","SVTYPE","SVLEN","TOOL","PAFV","PAFG","PAF1KG","GC","CR","MQ","SEGD","NUMG","GENES","GFEAT","HPO","PHEN");
@columnsMin=("SAMPLE","ID","FT","SU","DRF","DRA","CNV","IGV","GOTO","LOCATION","SVTYPE","SVLEN","TOOL","PAFV","PAFG","NUMG","GENES","GFEAT","PHEN");
# @columns=("SAMPLE","FT","RARE","SU","PE","SR","DRF","DRA","IDD","LOCATION","SVTYPE","SVLEN","TOOL","PAFV","PAFDRA","PAF1KG","GC","NUMG","GENES","HPO"); ## Cooper version

if($ped){
	undef @t;  map{ push @t,"pedInfo$_" } 1..$infoFields;  @columns=( "family", @t, "affected","IA","IUA",@columns ); @columnsMin=( "family", @t, "affected","IA","IUA",@columnsMin ); #"SUC","DRAC",
# 	$V->add_header("##FORMAT=<ID=SUC,Number=1,Type=Float,Description=\"Frequency of SR+PE relative to avg of unaffacted individuals. 0: absent in individual, 1: absent in unaffacted, 0.5: half as average unaffacted; -1 no value for unaffacted \">"); 
# 	$V->add_header("##FORMAT=<ID=DRAC,Number=1,Type=Float,Description=\"Frequency of DRA relative to avg of unaffacted individuals. 0: absent in individual, 1: absent in unaffacted, 0.5: half as average unaffacted; -1 no value for unaffacted \">"); 
	$V->add_header("##INFO=<ID=IA,Number=1,Type=Integer,Description=\"Number of times a variant was detected in affected individuals\">"); 
	$V->add_header("##INFO=<ID=IUA,Number=1,Type=Integer,Description=\"Number of times a variant detected in unaffected individuals\">"); 
}
if($geneL){ push @columns, "CANDG" ; push @columnsMin, "CANDG" ;}

print VCF ${$V->print_header};
print RAREPV ${$V->print_header};
print PASS ${$V->print_header};
$rowC=0;
while (my $V=$V->next_line){
	
	$c++; 
		
	next if exists($$V{INFO}{SECONDARY});
	
	if(exists($$V{INFO}{CNEUTR}) and !exists($$V{INFO}{CNV})){  $$V{INFO}{CNV}=($$V{INFO}{CNEUTR}==1)? 0:1;  }
	
	
	$p=0; $r=0; $passTrustedTest=0;
	foreach $cSample ( @{$$V{samples}} ){
	 	 
		### rare filter
		$$V{GTS}{RARE}{$cSample}=1;

		# in case of PKD1 first exon does shows up as rare in PCR free, even though SU PAFSU clearly show that it is rare. Better trust PAFSU if SU sufficiently high else PAFDRA
		if( exists($$V{GTS}{SU}{$cSample}) and exists($$V{GTS}{PAFSU}{$cSample}) ){
			$$V{GTS}{RARE}{$cSample}=0 if $$V{GTS}{PAFSU}{$cSample}>0.01;# and $$V{GTS}{PAFSU}{$cSample}<1.4);
			$passTrustedTest++;	 
		}
		if(exists($$V{GTS}{PAFDRA}{$cSample}) and $$V{INFO}{CNV}==1 and $$V{GTS}{PAFDRA}{$cSample}>=0 ){
			$$V{GTS}{RARE}{$cSample}=0 if $$V{GTS}{PAFDRA}{$cSample}>0.01;
			$passTrustedTest++;	  
		}

		if(exists($$V{INFO}{PAFV}) and $$V{INFO}{PAFV}>=0 ){
			$$V{GTS}{RARE}{$cSample}=0 if $$V{INFO}{PAFV}>0.01;
			$passTrustedTest++;	 
		}
		if(exists($$V{INFO}{PAFG}) and $$V{INFO}{PAFG}>=0 ){
			$$V{GTS}{RARE}{$cSample}=0 if $$V{INFO}{PAFG}>0.01;
			$passTrustedTest++;	 
		}
		if( $passTrustedTest == 0 and exists($$V{INFO}{PAF1KG}) and $$V{INFO}{PAF1KG}>=0 ){
			$$V{GTS}{RARE}{$cSample}=0 if $$V{INFO}{PAF1KG}>0.01;
		}


# 		$$V{GTS}{RARE}{$cSample}=2 if $$V{GTS}{RARE}{$cSample}<=0 and length(getVal("PHEN"))>0; ### include all CNVs overlapping pathogenic gene

		$$V{GTS}{RARE}{$cSample}=(-1) if $$V{GTS}{PAFDRA}{$cSample} == (-1) and $$V{GTS}{PAFSU}{$cSample} == (-1) and $$V{INFO}{PAF1KG} == (-1) and $$V{INFO}{PAFV} < 0 ;

		$r++ if $$V{GTS}{RARE}{$cSample}!=0;


		next if exists($$V{GTS}{FT}{$cSample}) and $$V{GTS}{FT}{$cSample} ne "PASS"; 
		$$V{GTS}{FT}{$cSample}="LOW";
		next if exists($$V{INFO}{LOCATION}) and $$V{INFO}{LOCATION} =~ /hs37d5/;

		next if $chrPf eq '' and $$V{CHROM} !~ /^(X|Y|MT|[0-9]+)$/ ;
		next if $chrPf eq 'chr' and $$V{CHROM} !~ /^(chrX|chrY|chrM|chr[0-9]+)$/ ;
		
		
		
		if($$V{INFO}{CNV} == 0 ){
		  next if $$V{GTS}{MQBP}{$cSample}<50 and $$V{GTS}{SU}{$cSample}<6;
		  next if $$V{INFO}{SVTYPE} eq "BND" and $$V{GTS}{SU}{$cSample}<6;

		}else{
		 if ($$V{GTS}{SU}{$cSample}<=1 and $$V{INFO}{SVLEN}<10000){
			next;
		 }
		 if($$V{INFO}{SVTYPE} eq "DEL"){
			next if exists($$V{GTS}{DRF}{$cSample}) and $$V{GTS}{DRF}{$cSample}>=0.8 and $$V{GTS}{DRA}{$cSample}>=0.8;
		 }elsif($$V{INFO}{SVTYPE} eq "DUP"){
			next if exists($$V{GTS}{DRF}{$cSample}) and $$V{GTS}{DRF}{$cSample} <= 1.2 and $$V{GTS}{DRA}{$cSample} <= 1.2; 
		 }	
		}	 
		
		
		
		#print STDERR $$V{GTS}{PAFDRA}{$cSample}." ".$$V{GTS}{PAFSU}{$cSample}." ".$$V{INFO}{PAF1KG}."\n" if $$V{GTS}{PAFDRA}{$cSample} == (-1) and $$V{GTS}{PAFSU}{$cSample} == (-1) and $$V{INFO}{PAF1KG} == (-1);

		$$V{GTS}{FT}{$cSample}="PASS";$p++;
		if($$V{GTS}{SU}{$cSample}>10 & $$V{GTS}{SR}{$cSample}>0 & $$V{GTS}{PE}{$cSample}>0){ $$V{GTS}{FT}{$cSample}="HIGH"; }
		if(($$V{INFO}{SVTYPE} eq "DEL" or $$V{INFO}{SVTYPE} eq "DUP") & $$V{INFO}{CNV} == 1 & ($$V{INFO}{SVLEN}>100000 or ($$V{INFO}{SVLEN}>10000 and exists($$V{INFO}{MQ}) and $$V{INFO}{MQ}>55 )) ){ $$V{GTS}{FT}{$cSample}="HIGH"; }
	 
	}
	
	addPedInfo() if $ped;
	if ($geneL){  $cS=$$V{INFO}{GENES}; $cS=~ s/"//g;  $$V{INFO}{CANDG}="\"".join(",",     grep { exists($geneList{$_}) } split(",",$cS) )."\"";   } # mark genes in gene list
	
	$$V{FILTER}="PASS" if $p>0;
	$$V{FILTER}="LOW" if $p==0;
	
	print VCF ${$V->print_line};
	print RAREPV ${$V->print_line} if $r>0 and $p>0 and $$V{INFO}{NUMG}>0;
	print PASS ${$V->print_line} if $p>0;
	
	#### output text
	if ( $c==1 ){ 
		if(exists($$V{GTS}{HR}{$cSample})){  push @columns,"HR" ; push @columns,"HC"; } 
		print TXT join("\t",@columns)."\n"; print RAREP join("\t",@columns)."\n";  
		$xl{ws}{"full"}->write_row($rowC, 0, \@columns);
		$rowC2++;
		$xl{ws}{"RGfull"}->write_row($rowC, 0, \@columns);
		$xl{ws}{"RGlight"}->write_row($rowC, 0, \@columnsMin);
		$rowC++;
	}
	
	foreach $cSample ( @{$$V{samples}} ){
		undef @cOut; map {  push @cOut, getVal($_);  } @columns;
		undef @cOutMin; map {  push @cOutMin, getVal($_);  } @columnsMin;
		print TXT join("\t",@cOut)."\n";
		
		$xl{ws}{"full"}->write_row($rowC2, 0, \@cOut);
		$rowC2++;
		
		if ( $$V{GTS}{RARE}{$cSample}!=0 and ($$V{GTS}{FT}{$cSample} eq "PASS" or $$V{GTS}{FT}{$cSample} eq "HIGH") and $$V{INFO}{NUMG}>0){
			print RAREP join("\t",@cOut)."\n";
			$xl{ws}{"RGfull"}->write_row($rowC, 0, \@cOut);
			$xl{ws}{"RGlight"}->write_row($rowC, 0, \@cOutMin);
			$rowC++;
		}


	}
	
}


# format excel worksheet
################

for (my $i=0; $i<=$#columnsMin; $i++){ $xl{ws}{"RGlight"}->set_column($i, $i, 6.5); }
for (("RGfull","full")){
	for (my $i=0; $i<=$#columns; $i++){ $xl{ws}{$_}->set_column($i, $i, 5.5); }
}

for (("RGfull","RGlight","full")){
	$xl{f1}{$_} = $xl{wb}{$_}->add_format();
	$xl{f1}{$_}->set_num_format('#,##0');
	$xl{f2}{$_} = $xl{wb}{$_}->add_format();
	$xl{f2}{$_}->set_num_format("\@");
}



%col2width=("SAMPLE" => 10.5, "SVTYPE" => 4.5, "PHEN" => 8,"TOOL" => 8.8, "SVLEN" => 10, "HPO" => 6.5);
%col2format=("SVLEN" => $xl{f1} ,"HPO" => $xl{f2});

for (("RGfull","full")){
	format_excl(\@columns, $xl{ws}{$_}, $_ );
}
format_excl(\@columnsMin, $xl{ws}{"RGlight"}, "RGlight" );


sub format_excl{

	my( $columns, $cws, $fullOrMin )=@_;
	my %columns2num;
	
	for (my $i=0; $i<=$#$columns; $i++){ $columns2num{$$columns[$i]}=$i; }
	
	foreach (keys %col2width){
		next if !exists($columns2num{$_});
		$cCol=$columns2num{$_};
		$cws->set_column($cCol, $cCol);
	}	
	foreach (keys %col2format){ # add format
		next if !exists($columns2num{$_});
		$cCol=$columns2num{$_};
		$cws->set_column($cCol, undef, $col2format{$_}->{$fullOrMin});
	}	
}

for (("RGfull","RGlight","full")){
	$xl{ws}{$_}->freeze_panes(1, 0);
	$xl{wb}{$_}->close();
}


close(TXT);
close(VCF);
close(RAREP);
close(RAREPV);
close(PASS);


sub getVal {

	($cCol)=@_;
	
	if($cCol eq "PHEN"){
		undef @tmpTxt;
		$cGeneList=$$V{INFO}{GENES};
		$cGeneList=~ s/"//g;
		foreach $cGene (split(",",$cGeneList)){
			if ( exists($ph_gene2phen{$cGene}) ){
				push @tmpTxt, "$cGene: ".$ph_gene2phen{$cGene};
			}
		}
		$retVal=(@tmpTxt)? join("|",@tmpTxt):"";
# 		print STDERR "$retVal\n" if @tmpTxt;
		return $retVal;
	}elsif($cCol eq "IGV"){
		
		return "=HYPERLINK(\"http://localhost:60151/load?file=$IGVSessionDir/$cSample.xml&merge=false\",\"IGV\")";
		
	}elsif($cCol eq "GOTO"){
		
		
		if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
			if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			($cChr2,$cPos2)=($1,$2);
			}
		}else{
			($cChr2,$cPos2)=($$V{CHROM},$$V{"INFO"}{"END"});
		}
		
		
		$flan2=1000;
		
		if( $$V{CHROM} eq $cChr2 and $$V{INFO}{CNV}==1 ){ $flan=int(($$V{INFO}{"END"}-$$V{POS})*1.5);  $min1=($$V{POS}-$flan); $min1=1 if $min1<=0;  $cReg=$$V{CHROM}.":".$min1."-".($$V{INFO}{"END"}+$flan)   }
		else{  $min1=($$V{POS}-$flan2); $min1=1 if $min1<=0; $min2=($cPos2-$flan2); $min2=0 if $min2<0;  $cReg=$$V{CHROM}.":".$min1."-".($$V{POS}+$flan2)."%20".$cChr2.":".$min2."-".($cPos2+$flan2) }
		
		return "=HYPERLINK(\"http://localhost:60151/goto?locus=$cReg\",\"GOTO\")";
		
	}elsif($cCol eq "HPO"){
		
		if(exists($$V{INFO}{$cCol}) and $$V{INFO}{$cCol} !~ /^"[|]+"$/){ return unquote($$V{INFO}{$cCol}) }
		else{ return "" }
	}elsif($cCol eq "GFEAT"){
		# return the most sever part of any of the affected genes  ("start_codon","stop_codon","CDS","exon","transcript")
		
		
		
		if(exists($$V{INFO}{"GDESC"})){
		
			%tmp_h=();
			foreach ( split("[|]",$$V{INFO}{"GDESC"}) ){
				@tmp=split(",", $_ );
				$tmp_h{$tmp[2]}++;
				
			}
			for ( ("start_codon","stop_codon","CDS","exon","transcript")){
				
				if (exists($tmp_h{$_})){
					if ($_ eq "transcript"){
						return "intron"
					}else{
						return $_
					}
				}
				
			}
		}
		return ""
	}

	
	if($cCol eq "GENES" or $cCol eq "CANDG"){  return unquote($$V{INFO}{$cCol}) if exists($$V{INFO}{$cCol}); }
	
	return $cSample if $cCol eq "SAMPLE";
	
	return $sa2info{$cSample}{$cCol} if exists($sa2info{$cSample}{$cCol});
	return $$V{GTS}{$cCol}{$cSample} if exists($$V{GTS}{$cCol}{$cSample});
	return $$V{INFO}{$cCol} if exists($$V{INFO}{$cCol});	
	return $$V{$cCol} if exists($$V{$cCol});	
	return ".";
}





sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}
sub roundA{ ($cFreq)=@_;  if($cFreq!=1){ $RL=int(log10($cFreq)*(-1))+2; $cFreq=round($cFreq,$RL);  }  return $cFreq }
sub log10 { my $n = shift;  my $r=($n>0)? log($n)/log(10):0; return $r; }


sub addPedInfo{

	
	$SUC=0;$SUCC=0; 
	$DRAC=0;$DRACC=0; 
	$$V{INFO}{IA}=0;
	$$V{INFO}{IUA}=0;

	
	foreach $cSample ( @{$$V{samples}} ){
		if (exists($sa2nonAff{$cSample})){
# 			print "$cSample,";
			if($$V{GTS}{DRA}{$cSample} ne "."){ $DRAC+=$$V{GTS}{DRA}{$cSample};  $DRACC++ }
			
			if($$V{GTS}{SU}{$cSample} ne "."){ $SUC+=$$V{GTS}{SU}{$cSample};  $SUCC++ }
			
			$$V{INFO}{IUA}++ if ($$V{GTS}{FT}{$cSample} eq "PASS" or $$V{GTS}{FT}{$cSample} eq "HIGH");
			
# 			print "$cSample: SU:$$V{GTS}{SU}{$cSample},  DRA:$$V{GTS}{DRA}{$cSample}\n";
		}
		if (exists($sa2Aff{$cSample})){
			$$V{INFO}{IA}++ if ($$V{GTS}{FT}{$cSample} eq "PASS" or $$V{GTS}{FT}{$cSample} eq "HIGH");
		}
	}
	
# 	foreach $cSample ( @{$$V{samples}} ){
# 		
# 		if($SUCC>0){
# 			$$V{GTS}{SUC}{$cSample}=($$V{GTS}{SU}{$cSample}>0)? roundA(($SUC/$$V{GTS}{SU}{$cSample})/$SUCC):(-1); # frequence relative to average onAffected
# # 			print "$cSample ($SUC / ".$$V{GTS}{SU}{$cSample}.") / $SUCC    $_\n";
# 		}else{
# 			$$V{GTS}{SUC}{$cSample}=(-1);  # can not be calculated because no non affected has a value 
# 		}
# 		
# 		if($DRACC>0){
# 			$$V{GTS}{DRAC}{$cSample}=($$V{GTS}{DRA}{$cSample}>0)? roundA(($DRAC/$$V{GTS}{DRA}{$cSample})/$DRACC):100; # frequence relative to average onAffected 
# 		}else{
# 			$$V{GTS}{DRAC}{$cSample}=(-1);  # can not be calculated because no non affected has a value
# 		}
# # 			print "    $cSample  SU:$$V{GTS}{SU}{$cSample}   SUC:$$V{GTS}{SUC}{$cSample},  DRA:$$V{GTS}{DRA}{$cSample}: DRAC:$$V{GTS}{DRAC}{$cSample}\n";		
#  		
#  	}


}

${$BrStat{$cBrID}{PAFSU}}[$cID]=($cSU>0)? roundA(  ($BrStat{$cBrID}{SUC}/$cSU)/$S_SV_control_numSamples   ):(-1);  






sub unquote{
	
	my ($s)=@_;
	if(substr($s,0,1) eq "\"" and substr($s,-1,1) eq "\"" ){
		return substr($s, 1, -1)
	}else{
		return $s;
	}

}	

































