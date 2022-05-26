
# Author: Andre Minoche, a.minoche@garvan.org.au
# 


package My::VCF2;
use strict;
use warnings;
 

my @header;
my @colNames=("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");


my %TrueFlaseInfoFields;

my %settings;

my $headLineC=0;

############# initiate new object
# my $IntersObj=My::InterS->new;
sub new{
    my $class = shift;
    my $inF= shift;
#     my $headLineC= 0;
    map {$settings{$_}++} @ARGV; 
    
    my $s = {
        header => {},
        headerOrder => {}
    };
    bless $s, $class;
    
	my $cHeaderType="";
	if($inF =~ /[.]gz$/){open($$s{fhandle}, "gzip -dc $inF |") || die "can not open $inF";}
	else{open($$s{fhandle}, "<$inF") || die "can not open $inF";}
	
	while(readline($$s{fhandle})){
		chomp; 
		die "missing header $_" if $_ !~ /^#/;
		$headLineC++;


				
		if(/^##([^=]+)=/ or /^#(CHROM)/){
			$cHeaderType=$1;
			
			$$s{headerOrder}{$cHeaderType}=($cHeaderType eq "CHROM") ? 1000000:$headLineC;
			push @{$$s{header}{$cHeaderType}}, $_;
# 			print "new $cHeaderType ".$$s{headerOrder}{$cHeaderType}." $_\n";
			if (/^##INFO=<ID=([^,]+),Number=0/){ $TrueFlaseInfoFields{$1}++; }
			
			if (/^##(INFO)=<ID=([^,]+)/){ 
					$$s{TagOrderC}{$1}++; $$s{TagOrder}{$1}{$2}=$$s{TagOrderC}{$1}; 
			}
# 			if (/^##(INFO|FORMAT)=<ID=([^,]+)/){ $$s{TagOrderC}{$1}++; $$s{TagOrder}{$1}{$2}=$$s{TagOrderC}{$1}; }
			
		}else{
			die "odd header 2 $_\n";
		}
		
		
		
		if(/^#CHROM/){
			my @L=split("\t",$_);
			# print STDERR "$_\n"; exit;
			for (my $j=9; $j<=$#L;$j++){ push @{$$s{samples}}, $L[$j]; $$s{sample2pos}{$L[$j]}=$j; $$s{pos2sample}{$j}=$L[$j];     }
			last;
		}		
		

		
	}	
    return $s;
}

sub remove_header{

	my $s = shift;
    my $headerTypeToDelete = shift;
	my $headerToDelete = shift;
	
	$$s{skip}{$headerTypeToDelete}{$headerToDelete}++;

}


sub get_header_type{
	
	my $header=shift;
	
	my ($cHeaderType, $headerID) = ("","");
	
	if($header =~ /^##([^=]+)=/){
		$cHeaderType=$1;
	}else{die "odd header 3 $header\n";}
	
	# if the header contains an ID, check if ID already exists, if so, die
	if($header =~ /^##[^=]+=<ID=([^,]+),/){
		$headerID=$1;
		
	}elsif($header =~ /^##[^=]+=[^,]+$/){
		#print STDERR " this is a short header: $cHeaderType $header\n";		
		$headerID=$cHeaderType;
	}else{
		die "odd header 4 $header\n";
	}
	
	return ($cHeaderType, $headerID)
	
}

sub add_header{
    my $s = shift;
    my $headerToAdd = shift;

	my ($cHeaderType,$headerToAddID) = get_header_type($headerToAdd);
	
	# overwrite header if already present
	my $header_replaced=0;
	my ($cHeaderType2,$headerToAddID2) = ("","");
	for(my $i=0; $i<=$#{$$s{header}{$cHeaderType}}; $i++){
		
		($cHeaderType2,$headerToAddID2) = get_header_type(${$$s{header}{$cHeaderType}}[$i]);			
		if ($cHeaderType eq $cHeaderType2 and $headerToAddID eq $headerToAddID2){
			#print STDERR "\n - overwriting header from ".${$$s{header}{$cHeaderType}}[$i]." to ".$headerToAdd;
			${$$s{header}{$cHeaderType}}[$i]=$headerToAdd;
			$header_replaced=1;
		}	

	}
	
	if ($header_replaced==0){
		push @{$$s{header}{$cHeaderType}}, $headerToAdd; 

	
		if ( !exists($$s{headerOrder}{$cHeaderType} )){
			$$s{headerOrder}{$cHeaderType}=$headLineC;
			$headLineC+=1;
		}

		if ($headerToAdd =~ /^##(INFO|FORMAT|FILTER)=<ID=([^,]+)/){ 
			if( !exists($$s{TagOrder}{$1}{$2}) ){
				$$s{TagOrderC}{$1}++; $$s{TagOrder}{$1}{$2}=$$s{TagOrderC}{$1}; 
			}
		}elsif ($headerToAdd =~ /^##[^=]+=[^=]+$/){
			# no tag order needed to be recorded
		}else{
			die "odd header 5 $headerToAdd\n";
		}
	
	}
}

sub next_line{

	my $s = shift;
		
	while(readline($$s{fhandle})){
		
		$s=parse_line($s,\$_);
		
		### for freebayes
# 		$s=mostLikelyGT($s) if $settings{"mostLikelyGT"}; # correct freebayes GT bug 0/1:29:0:0:29:1065:-96.2172,-8.72987,0 (GT:DP:RO:QR:AO:QA:GL)
# 		$s=checkQualRefAlt($s) if $settings{"checkQualRefAlt"}; #qualRefAlt filter: for variant calls keep: "averge qual ref" - "average qual alt" < 10, if multiple alt present take the alt allele from GT
# 		$s=checkStrand($s) if $settings{"checkStrand"}; # strand fitler: for each VCF line: at least 10% of alt observations on reverse/forward strand, 
	
		### add
# 		$s=addAltCount($s) if $settings{"addAltCount"}; # count samples with alt allele >0 per vcf line
	
		return $s;
	}
	close($$s{fhandle});
	return 0;
}


sub parse_line {
		
		my ($s,$sLine)=@_;
		
		
		chomp($$sLine); 
		my @cLine=split("\t",$$sLine);
	    $$s{"all"}=$_;
	    
		foreach my $i ((0,1,2,3,4,5,6)){ $$s{$colNames[$i]}=$cLine[$i] } # CHROM-FILTER
		$$s{"CHROM-FILTER"}=join("\t",@cLine[0..6]);
		$$s{"INFO_all"}=$cLine[7];
		$$s{"FORMAT_all"}=$cLine[8];
		$$s{"GT_all"}=$cLine[9];
		
# 		print STDERR "fffffffff0   ".$#{$$s{samples}}."\n";
		
		for(my $i=0; $i<=$#{$$s{samples}}; $i++){
			my $str="GT_".${$$s{samples}}[$i];
# 			print "$str\n";
			$$s{$str}=$cLine[9+$i];
		}
		# parse INFO field
		undef %{$$s{INFO}};
		foreach my $p ( (split(";",$cLine[7])) ){
			my @x=split("=",$p); 
			$$s{INFO}{$x[0]}=(scalar(@x)>1)? $x[1]:1;
# 			print "$x[0] ".$$s{INFO}{$x[0]}."\n";
		}
		
		
		
		
		# parse GT field
		undef %{$$s{GTS}};
		my @x1=split(":",$cLine[8]); 
		for (my $j=9; $j<=$#cLine;$j++){  # for each sample
			my @x2=split(":",$cLine[$j]); 
			for (my $i=0; $i<=$#x1;$i++){  # for each GT field
				
				
				$$s{GTS}{$x1[$i]}{$$s{pos2sample}{$j}}=(exists($x2[$i]))? $x2[$i]:".";
				
				if( !exists($$s{TagOrder}{"FORMAT"}{$x1[$i]}) ){  $$s{TagOrderC}{"FORMAT"}++; $$s{TagOrder}{"FORMAT"}{$x1[$i]}=($$s{TagOrderC}{"FORMAT"}-1000); }
			} 
		}
		
		my $isIndel=0; map { $isIndel++ if length($_)>1} (split(",",$$s{ALT}),$$s{REF});
		$$s{is}{INDEL}=($isIndel>0)? 1:0;
		return $s;
		
}

########## reconstruct the VCF file
sub print_header{

	my ($s)=@_;
	my $cHeader="";
	
	
# 	print " order: ".join(' ',keys %{$$s{headerOrder}} )."\n";
# 	print " header: ".join(' ',keys %{$$s{header}} )."\n";
	
	foreach my $cHeaderType ( sort {$$s{headerOrder}{$a} <=> $$s{headerOrder}{$b} } keys %{$$s{headerOrder}} ){
		
		
		foreach my $cHeaderLine (@{$$s{header}{$cHeaderType}}){
		    
		    #print STDERR "print_header $cHeaderType ".$$s{headerOrder}{$cHeaderType}." $cHeaderLine\n";
		    
		    if($cHeaderLine =~ /^##[^=]+=<ID=([^,]+)/){ next if exists($$s{skip}{$cHeaderType}{$1}); }
			
			$cHeader.=$cHeaderLine."\n";
		}
	}
	return \$cHeader;		
}


sub print_line{

	my ($s)=@_;
	my @cLine=();
	my $outLine="";
	foreach my $i ((0,1,2,3,4,5,6)){ $cLine[$i]=$$s{$colNames[$i]}; } # CHROM-FILTER
	
	# INFO
	$cLine[7]="";
	foreach my $cI ( sort { $$s{TagOrder}{INFO}{$a} <=> $$s{TagOrder}{INFO}{$b} } keys %{$$s{TagOrder}{INFO}} ){ 
		next if !exists($$s{INFO}{$cI});
		next if exists($$s{skip}{INFO}{$cI});
		
		if(exists($TrueFlaseInfoFields{$cI})){ 
			$cLine[7].=(exists($$s{INFO}{$cI}))? "$cI;":"" 
		}else{ 
			if (!exists($$s{INFO}) or !exists($$s{INFO}{$cI}) or !defined $cI or !defined $$s{INFO}{$cI}){ print STDERR " ! field $cI does not exist in ".$$s{all}."\n"; }
			$cLine[7].="$cI=".$$s{INFO}{$cI}.";" 	 
		}
		
	} chop($cLine[7]); 
	
	# FORMAT
	foreach my $cI ( sort { $$s{TagOrder}{FORMAT}{$a} <=> $$s{TagOrder}{FORMAT}{$b} } keys %{$$s{TagOrder}{FORMAT}} ){ 
		next if !exists($$s{GTS}{$cI});
		next if exists($$s{skip}{FORMAT}{$cI});
		$cLine[8].="$cI:";
	} chop($cLine[8]); 
		
	# GT	
	foreach my $cSample (keys %{$$s{sample2pos}}){ # foreach sample
		my $cPos=$$s{sample2pos}{$cSample};
		foreach my $cI ( sort { $$s{TagOrder}{FORMAT}{$a} <=> $$s{TagOrder}{FORMAT}{$b} } keys %{$$s{TagOrder}{FORMAT}} ){ # foreach GT field
			next if !exists($$s{GTS}{$cI}); # current lines does not contain all GT fields
			next if exists($$s{skip}{FORMAT}{$cI});
			$cLine[$cPos].=(exists($$s{GTS}{$cI}{$cSample}))? $$s{GTS}{$cI}{$cSample}.":":"0:";  
		}
		chop($cLine[$cPos]);
	}
	$outLine=join("\t",(@cLine))."\n";
	return \$outLine;
		
}



1;
