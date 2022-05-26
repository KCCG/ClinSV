
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 

($inGenes,$inVCF,$S_HPO_gene2label,$S_HPO_gene2hpo)=@ARGV;

%allowedChr;

$verbose=0;
$maxGenesPerVariant=100;




############################# create hash with gene symbol to HPO
%acc2symb;
%geneSymb2HPO;

open(IN2, "<$S_HPO_gene2label") || die "$S_HPO_gene2label nicht gefunden";
while(<IN2>){ chomp; @_=split(/[\t:]/,$_);  $acc2symb{$_[1]}=$_[2]; #print "$_[1] $_[2]\n";
}close(IN2);

open(IN1, "<$S_HPO_gene2hpo") || die "$S_HPO_gene2hpo nicht gefunden";
while(<IN1>){ chomp; @_=split(/[\t:]/,$_); push @{$geneSymb2HPO{$acc2symb{$_[1]}}}, $_[3] if exists($acc2symb{$_[1]}); # print $acc2symb{$_[1]}." $_[3]\n";
}close(IN1);

print STDERR scalar(keys %geneSymb2HPO)." genes with HPO annotation\n";

############################# annotate the VCF file







print STDERR "parsing gff...\n";
##################################################################
open(IN1, "zcat $inGenes |") || die " nicht gefunden";
while(<IN1>){ chomp; @_=split("\t",$_); 
	$allowedChr{$_[0]}++; 
	$c++;
}close(IN1);

close(IN1);
$NrAllSamples=scalar(keys %samplesH);
print STDERR "input features $c on ".scalar(keys %allowedChr)." reference sequences\n";



print STDERR "annotating vcf...\n";
##################################################################
use My::VCF2; 
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
my $tabix = new Tabix(-data => "$inGenes",-index => "$inGenes.tbi");


my $VCFObj=My::VCF2->new($inVCF);
$VCFObj->add_header("##INFO=<ID=NUMG,Number=1,Type=Integer,Description=\"Number of genes affected by the SV\">");
$VCFObj->add_header("##INFO=<ID=GENES,Number=1,Type=String,Description=\"ENSEMBL genes affected by the variant\">");
$VCFObj->add_header("##INFO=<ID=GDESC,Number=1,Type=String,Description=\"Details about affected gene: gene name, gene type, affected feature, transcript ID\">");
$VCFObj->add_header("##INFO=<ID=HPO,Number=1,Type=String,Description=\"HPO numbers of affected genes. HPO's of genes are separated by | and appear in the same order as the gene names in the GENES column. Multiple HPO's per gene are separated by colons. \">");
print ${$VCFObj->print_header};

$variantsAffectedByGene=0;
$allVariantsInVCF=0;

while (my $V=$VCFObj->next_line){
  		
  	$allVariantsInVCF++;
  	undef %affectedGenesN;
	undef @affectedGenesD;
	undef %GTF;	
	
  	($cChr1,$cPos1)=($$V{CHROM},$$V{POS});
  	
  	print ${$V->print_line} if $verbose;
	
	if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
# 		next if exists($$V{"INFO"}{"SECONDARY"});
		if($$V{"ALT"} =~ /[\[\]](.+):([0-9]+)[\[\]]N*$/){
			
			($cChr2,$cPos2)=($1,$2);
			if($cChr1 eq $cChr2){ 	
				$$V{"INFO"}{"LEN"}=abs($cPos1-$cPos2);
				$$V{"INFO"}{"END"}=$cPos2;
			}
		}else{die $$V{"all"} }	
	}else{
		($cChr2,$cPos2)=($$V{CHROM},$$V{"INFO"}{"END"});
	}
		

	######### if copy neutral only check at breakpoint, else annotate everything in between	
	if(  (exists($$V{"INFO"}{"CNV"}) and $$V{"INFO"}{"CNV"}==1) ){ # or  $$V{"INFO"}{"SVTYPE"} eq "INV" ){
		
		######### check region overlap
		checkReg($cChr1,$cPos1,$cChr2,$cPos2);
		
	}else{
		
		######### check breakpoint
		checkBR($cChr1,$cPos1,$cChr2,$cPos2);
		
	}
	
	$affectedGenes=scalar(keys %affectedGenesN);
	$$V{"INFO"}{"NUMG"}=$affectedGenes;
	if($affectedGenes>0 and $affectedGenes<=$maxGenesPerVariant){
		$variantsAffectedByGene++;
		$$V{"INFO"}{"GENES"}="\"".join(",",(sort keys %affectedGenesN))."\"";
		
		# genereate HPO string
		undef (@affectedHPOs);
		foreach $cSymb (sort keys %affectedGenesN){
			$cHPOs=(exists($geneSymb2HPO{$cSymb}))? $geneSymb2HPO{$cSymb}:"";
			push @affectedHPOs, join(":",@$cHPOs);
		}
		$$V{"INFO"}{"HPO"}="\"".join("|",(@affectedHPOs))."\"";
		
		reportMostImportantFeature();
		$$V{"INFO"}{"GDESC"}="\"".join("|",@affectedGenesD)."\"";
	}elsif($affectedGenes>$maxGenesPerVariant){
		$variantsAffectedByGene++;
		$$V{"INFO"}{"GENES"}="\"more than $maxGenesPerVariant genes...\"";
		reportMostImportantFeature();
		$$V{"INFO"}{"GDESC"}="\"more than $maxGenesPerVariant genes...\"";
	}
	
	print ${$V->print_line};
	
	# 	last if $c2++>100;
	
}

sub checkBR{
		
		my ($cChr1,$cPos1,$cChr2,$cPos2)=@_;
		@testCoords=([$cChr1,$cPos1,$cChr2,$cPos2],[$cChr2,$cPos2,$cChr1,$cPos1]);
		
		foreach $cA (@testCoords){
			
			my ($cChr,$cPos,$oChr,$oPos)=@$cA;
			
			if(!exists($allowedChr{$cChr}) ){return 0}
			print "--checkBR $cChr,$cPos,$oChr,$oPos\n" if $verbose;
		
			$res = $tabix->query($cChr,($cPos-5),($cPos+5)); 
			# 1. if only one breakpoint hits a gene, then it is affected by SV
			# 2. if both hit the gene but not the same intron, then the gene is also affected
			# in other words (since introns are not in gff), if overlap with any feature but the gene or transcript feature, then its affected
			if (defined $res->get){
			while(my $line = $tabix->read($res)){ #1       lincRNA exon    29554   30039   .       +       .       gene_id="ENSG00000243485";transcript_id="ENST00000473358";exon_number="1";gene_name="MIR1302-10";gene_source="ensembl_havana";gene_biotype="lincRNA";transcript_name="MIR1302-10-001";transcript_source="havana";exon_id="ENSE00001947070";
				
				chomp($line); ($bChr,$bSrc,$bFType,$bSt,$bEn,$bX1,$bX2,$bX3,$bName)=split("\t",$line);
				undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=substr($cTMP2[1],1,(length($cTMP2[1])-2)); }
				
				
				next if $bFType eq "gene" or $bFType eq "Selenocysteine";
				die "!transcript_id $line\n" if !exists($h{transcript_id});
				die "!gene_id $line\n" if !exists($h{gene_id});
				die "!gene_name $line\n" if !exists($h{gene_name});
				
				# remaining: CDS,UTR,exon,start_codon,stop_codon
				print "details hit BR:  ".scalar(keys %GTF)." $bFType ".$h{gene_name}." ".$h{transcript_id}." St,En:$bSt,$bEn\n" if $verbose;
				$affectedGenesN{$h{gene_name}}++;
				$GTF{$h{gene_id}}{$h{transcript_id}}{$bFType}=$h{gene_name}.",$bSrc,$bFType";
				return if scalar(keys %GTF)>$maxGenesPerVariant;
					
			}}
		}
		
}


sub checkReg{
		
		my ($cChr1,$cPos1,$cChr2,$cPos2)=@_;
		
		($cPos1,$cPos2)=($cPos1<$cPos2)? ($cPos1,$cPos2):($cPos2,$cPos1);
		print "--checkReg $cChr1,$cPos1,$cChr2,$cPos2\n" if $verbose;
		
		my ($qChr,$qPos1,$qPos2)=($cChr1,($cPos1-5),($cPos2+5));
		
		if(!exists($allowedChr{$qChr}) ){return 0}
		
		$res = $tabix->query($qChr,$qPos1,$qPos2); 
	
		# 1. if any of these features "CDS,UTR,exon,start_codon,stop_codon" are hit, the gene is affected
		while(my $line = $tabix->read($res)){ #1       lincRNA exon    29554   30039   .       +       .       gene_id="ENSG00000243485";transcript_id="ENST00000473358";exon_number="1";gene_name="MIR1302-10";gene_source="ensembl_havana";gene_biotype="lincRNA";transcript_name="MIR1302-10-001";transcript_source="havana";exon_id="ENSE00001947070";
			
			chomp($line); ($bChr,$bSrc,$bFType,$bSt,$bEn,$bX1,$bX2,$bX3,$bName)=split("\t",$line);
			undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=substr($cTMP2[1],1,(length($cTMP2[1])-2)); }
			
			next if $bFType eq "gene" or $bFType eq "Selenocysteine"; # $bFType eq "transcript" or 
			die "!transcript_id $line\n" if !exists($h{transcript_id});
			die "!gene_id $line\n" if !exists($h{gene_id});
			die "!gene_name $line\n" if !exists($h{gene_name});
					
			# remaining: CDS,UTR,exon,gene,start_codon,stop_codon
			print "details hit REG:  ".scalar(keys %GTF)." $bFType ".$h{gene_name}." ".$h{transcript_id}."\n" if $verbose;
			$GTF{$h{gene_id}}{$h{transcript_id}}{$bFType}=$h{gene_name}.",$bSrc,$bFType";
			$affectedGenesN{$h{gene_name}}++;
			
			return if scalar(keys %GTF)>$maxGenesPerVariant;
		}
}

	
sub reportMostImportantFeature{
	
	print "reportMostImportantFeature:  ".scalar(keys %GTF)."\n" if $verbose;
	#### report the most important feature that was hit		
	foreach $cGID (sort keys %GTF){
		
		$thisGeneAffected=0;
		foreach $cTID (sort keys %{$GTF{$cGID}}){
			
			# determine the most important feature
			# @FTypeImportanceOrder=("start_codon","CDS","stop_codon","exon","UTR","gene");
			foreach $cFType (("start_codon","CDS","stop_codon","exon","UTR","transcript","gene")){
			
				if(exists($GTF{$cGID}{$cTID}{$cFType})){
					print "reportMostImportantFeature2:  ".$GTF{$cGID}{$cTID}{$cFType}.",$cTID\n" if $verbose;
					# gene name, affected feature, TID
					push @affectedGenesD, $GTF{$cGID}{$cTID}{$cFType}.",$cTID";
					$thisGeneAffected=1;
					last;
				}
			}
			last if $thisGeneAffected; # report only one transcript per gene
		}
	}
}

	

print STDERR "of all variants in VCF $allVariantsInVCF, $variantsAffectedByGene (".sprintf("%.1f",($variantsAffectedByGene/$allVariantsInVCF*100))."%) hit ENSEMBL genes\n";



sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}







