
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";


use My::InterS;
use My::VCF2; 
use List::BinarySearch::XS qw( binsearch_pos );

($inSVs,$inSNVs)=@ARGV;


############################# annotate the VCF file

%allowedGT=("1/1",1,"0/1",1);

my $V=My::VCF2->new($inSVs);

$V->add_header("##FORMAT=<ID=HR,Number=1,Type=Float,Description=\"Number heterozygous variants divided by homozygous and heterozygous variants. HET:= VAF >0.11 and <0.89 \">");
$V->add_header("##FORMAT=<ID=HC,Number=1,Type=Integer,Description=\"Count heterozygous and homozygous variants within SV/CNV \">");
print ${$V->print_header};

$W=My::VCF2->new($inSNVs);


%h;

print STDERR " annotete SVs ... \n";

while (my $V=$V->next_line){
  	
  	 ($cChr1,$cPos1)=($$V{CHROM},$$V{POS});
	
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
	
	($cPos1,$cPos2)=sort {$a <=> $b} ($cPos1,$cPos2);
# 	next if $cChr1 ne "22";
	
	if ($cChr1 ne $cChr2){print ${$V->print_line}; next}
#   	if(!exists($$V{"INFO"}{"SVTORIG"}) or  # only check the copy neutral events
#   		($$V{"INFO"}{"SVTORIG"} ne "DUP" and $$V{"INFO"}{"SVTORIG"} ne "DEL" )  ){print ${$V->print_line}; next}
  	
#   	print STDERR "----- $cChr1,$cPos1,$cPos2 $$V{INFO}{SVLEN} \n";
  	
  	# read coords of current chromosome
  	if($chrCheck{$cChr1}++ == 0 ){ undef %h; readSNVsForChrom(); }
  	 	
  	foreach $cSample (@{$$V{samples}}){
  		
  		# count HET vars
		($iSt,$iEn)=binsearch_range($h{$cSample}{HET},$cPos1,$cPos2); 
		if($iSt>=0){
# 			print STDERR "->$iSt,$iEn $cChr1,$cPos1,$cPos2 ".${$h{$cSample}{HET}}[$iSt]."-".${$h{$cSample}{HET}}[$iEn]."\n"; 
			$cHET=$iEn-$iSt+1;
		}else{ $cHET=0; }
		
		
		# count HOM vars
		($iSt,$iEn)=binsearch_range($h{$cSample}{HOM},$cPos1,$cPos2); 
		if($iSt>=0){ 
# 			print STDERR "->$iSt,$iEn $cChr1,$cPos1,$cPos2 ".${$h{$cSample}{HOM}}[$iSt]."-".${$h{$cSample}{HOM}}[$iEn]."\n"; 
			$cHOM=$iEn-$iSt+1;
		}else{ $cHOM=0; }
		
		
  		# next if no variants called
  		if ($cHET==0 and  $cHOM==0){
  			$$V{GTS}{HC}{$cSample}=0;
  			next;
  		}
  		
  		$$V{GTS}{HR}{$cSample}=sprintf("%.0f",($cHET/($cHET+$cHOM)*100));
  		$$V{GTS}{HC}{$cSample}=($cHET+$cHOM);
  		
#   		print "W:  $cSample ".sprintf("%.0f",($cHET/($cHET+$cHOM)*100))." ".($cHET+$cHOM)."\n";
  	}

	print ${$V->print_line};
	
	
# 	last if $c2++>5;
}




sub readSNVsForChrom {

	print STDERR " read SNVs chrom $cChr1... \n";
	open(INVCF, "bcftools view -H -r $cChr1 $inSNVs  |") || die "can not open bcftools view $inSNVs $cChr1";
	while (<INVCF>){
	
		$W=$W->parse_line(\$_);
		
		next if $$W{FILTER} ne "PASS";
	
		foreach $cSample (@{$$W{samples}}){
				
			next if !exists( $allowedGT{$$W{GTS}{GT}{$cSample}} );
		
	# 			print "W:  $cSample $cBAF HET $h{$cSample}{HET} ".${$W->print_line}." \n";	
		
			@t=split(",",$$W{GTS}{AD}{$cSample});
			next if scalar(@t)!=2 or ($t[0]+$t[1])==0;
		
	# 			die "W:  $cSample $cBAF HET $h{$cSample}{HET} ".${$W->print_line}." \n" if ($t[0]+$t[1])==0;	
		
			$cBAF=($t[1]/($t[0]+$t[1]));

			if(0.11<=$cBAF and $cBAF<=0.89){
				push @{$h{$cSample}{HET}},$$W{POS};
	#   				print "W:  $cSample $cBAF HET $h{$cSample}{HET}\n";	
			
			}else{
				push @{$h{$cSample}{HOM}},$$W{POS};
	#   				print "W:  $cSample $cBAF HOM $h{$cSample}{HOM}\n";	
			}	
		
		}
	}
	
# 	print STDERR " sort coordinates ... \n";
	foreach $cSample (@{$$W{samples}}){
		@{$h{$cSample}{HET}}=sort {$a <=> $b} @{$h{$cSample}{HET}};
# 		print STDERR " HET $cSample ".scalar(@{$h{$cSample}{HET}})." coordinates ... \n";
		@{$h{$cSample}{HOM}}=sort {$a <=> $b} @{$h{$cSample}{HOM}};
# 		print STDERR " HOM $cSample ".scalar(@{$h{$cSample}{HOM}})." coordinates ... \n";
	}
	
	
	close(INVCF);
}


sub binsearch_range{ # examples 201,600->300-600; 200,200->200; 150,160->-1; 40-99->-1, 1001,1001->-1; 50,100->100, 1000,1000->1000, 1000,1001->1000, get all st and with same coord
 my ($r,$qSt,$qEn)=@_; ($qSt,$qEn)=sort {$a <=> $b} ($qSt,$qEn);
 $qSt=1 if $qSt<0; 
 if(($qSt<$$r[0] and $qEn<$$r[0]) or ($qSt>$$r[$#$r] and $qEn>$$r[$#$r])){ return ((((-1),(-1)))); } # St and En outside array
 my $iiSt = binsearch_pos {$a <=> $b} $qSt, @$r;
 my $iiEn = binsearch_pos {$a <=> $b} $qEn, @$r;
 $iiEn=$#$r if $iiEn>$#$r; $iiEn-- if $qEn<$$r[$iiEn] and $iiEn>0;
 while(($iiEn+1)<=$#$r and $$r[$iiEn]==$$r[$iiEn+1]){ $iiEn++; } # correct for multiple same entries upper end
 while(($iiSt-1)>=0 and $$r[$iiSt]==$$r[$iiSt-1]){ $iiSt--; } # correct for multiple same entries lower end
 if($iiSt>$iiEn){ return ((((-1),(-1))));} # 140,150 no interval in between
 return ($iiSt,$iiEn);
}






