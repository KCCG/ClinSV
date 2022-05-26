
# Author: Andre Minoche, a.minoche@garvan.org.au

# perl /g/data2/gd7/software/andmin/scripts/varpipe/sv/scripts/gene-centric-output.pl /g/data2/gd7/software/andmin/scripts/clinsv/refdata-b37/annotation/Homo_sapiens.GRCh37.75 /g/data2/gd7/research/andmin/tx70/R_160805_EMIMOU_LIONSDNA_M002/SVs/joined/tmp/SV-CNV.DGV.GC.LOH.SEGD.KDB.1kG.ENS.qsCNV.prioritized.vcf /g/data2/gd7/research/andmin/tx70/R_160805_EMIMOU_LIONSDNA_M002 $ref2 /g/data2/gd7/software/andmin/scripts/clinsv/refdata-b37/annotation/cosmic_gene_ID_list_Mar2017.ids > /g/data2/gd7/research/andmin/tx70/R_160805_EMIMOU_LIONSDNA_M002/SVs/joined/tmp/SV-CNV.DGV.GC.LOH.SEGD.KDB.1kG.ENS.qsCNV.prioritized.geneCentric.txt


use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2; 
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/^Subroutine Tabix.* redefined/)}; };
use Tabix;
use POSIX;
use Bio::DB::Big;
use File::Basename;

$verbose=1;


### for each variant, for each affected gene:

# gene name, gene ID, CANDG, gene type (coding, noncoding, miRNA) transcript ID, biotype, affected feature type (breakpoint location CDS, exon, intron) 
# affected exons exonA-exonB (covered by CNV, truncated) % of exon bases
# some SV fields IA, IUA, IGV,GOTO,  
# gene DRA (average of all exons, stdev)
# fusion partner (if linking 2 different genes):  geneB, fusion product: GeneA 5or3p exonA-exonB X% of exons pos, GeneB 5or3p exonA-exonB Y%, in frame, cosmic known fusion partner

$inGenes="/g/data2/gd7/research/andmin/resources/genes/Homo_sapiens.GRCh37.75";
($inGenes,$inVCF,$projectDir,$inRef,$geneListFile)=@ARGV;

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


$V=My::VCF2->new($inVCF);

# @annot_vcf_info=("SAMPLE","ID","FT","RARE","SU","PAFSU","PE","SR","DRF","DRA","DRAN","DRAP","GT","MQBP","CNV","IGV","GOTO","LOCATION","SVTYPE","SVLEN","TOOL","KDBF","VAF1KG","GC","CR","MQ","SEGD","NUMG","GENES","HPO","PHEN");

@annot_vcf_info=("ID","IA","IUA","SVLEN","TOOL","KDBF","CNV","IGV","GOTO","LOCATION","SVTYPE","SVLEN","GC","MQ");
@annot_vcf_gts=("SAMPLE","FT","RARE","SU","PE","SR","DRF","DRA","GT");
@outCatGene=("gene_name","gene_id","CANDG","biotype","tr_id","exonic_len","affected_feature");
@eCatGene=("exons_affected","exonic_bases_affected_pc","cds_bases_affected_pc","exons_affected_txt");
print join("\t",(  @outCatGene,@eCatGene,"fusion_gene_name",@eCatGene,@annot_vcf_info, (map { $_."_1" } @annot_vcf_gts ), (map { $_."_2" } @annot_vcf_gts )  ))."\n";





print STDERR "parsing gene list if available...\n";
##################################################################
if(-e $geneListFile){
	$geneL=1;
	print STDERR "    reading gene list $geneListFile ...\n";
	open(GLIST, "<$geneListFile") || die "can not write to $projectDir/sample.ped\n";
}elsif(-e "$projectDir/testGene.ids"){
	$geneL=1;
	print STDERR "    reading gene list $projectDir/testGene.ids ...\n";
	open(GLIST, "<$projectDir/testGene.ids") || die "can not write to $projectDir/sample.ped\n";
}
if($geneL){
	while(<GLIST>){		
		chomp; @t=split("\t",$_);
		$geneList{$t[0]}++;
# 			print "$t[0]\n";
	}
close(GLIST);
print STDERR "    number genes found:".scalar(keys %geneList)."\n";
}


print STDERR "parsing gene gff...\n";
##################################################################
%allowedChr;
open(IN1, "zcat $inGenes.gff.gz |") || die " nicht gefunden";
while(<IN1>){ chomp; @_=split("\t",$_); 
	$allowedChr{$_[0]}++; 
	$c++;
	
	($bChr,$bSrc,$bFType,$bSt,$bEn,$bX1,$bOri,$bX3,$bName)=split("\t",$_);
	undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=substr($cTMP2[1],1,(length($cTMP2[1])-2)); }
	
	# 1       lincRNA exon    30976   31097   .       +       .       gene_id="ENSG00000243485";transcript_id="ENST00000473358";exon_number="3";gene_name="MIR1302-10";gene_source="ensembl_havana";gene_biotype="lincRNA";transcript_name="MIR1302-10-001";transcript_source="havana";exon_id="ENSE00001827679";
	$bOri1=($bOri eq "+")? 1:(-1);
	if($bFType eq "gene"){
		$gid2{src}{$h{gene_id}}=$bSrc;
		$gid2{gene_name}{$h{gene_id}}=$h{gene_name};
		$gid2{gene_source}{$h{gene_id}}=$h{gene_source};
		$gid2{gene_biotype}{$h{gene_id}}=$h{gene_biotype};
	}elsif($bFType eq "transcript"){
		$tid2{orentation}{$h{transcript_id}}=$bOri1;
		$tid2{'chr'}{$h{transcript_id}}=$bChr;
		($tid2{transcript_st}{$h{transcript_id}}, $tid2{transcript_en}{$h{transcript_id}})=($bOri eq "+")? ($bSt,$bEn):($bEn,$bSt);
	}elsif($bFType eq "exon"){
		$tid2{exon_count}{$h{transcript_id}}++;
		$tid2{exonic_len}{$h{transcript_id}}+=($bEn-$bSt+1);
		$tid2{exons}{$h{transcript_id}}{$h{exon_number}}=[$bChr,$bSt,$bEn];	
# 		print STDERR  $h{transcript_id}." ".$bFType." ".$h{exon_number}." $bChr,$bSt,$bEn"." " .join("|",@{$tid2{exons}{ENST00000291386}{1}})."\n" if $h{transcript_id} eq "ENST00000291386";
	}elsif($bFType eq "CDS"){
		$tid2{cds_count}{$h{transcript_id}}++;
		$tid2{cds_len}{$h{transcript_id}}+=($bEn-$bSt+1);
		$tid2{cds}{$h{transcript_id}}{$h{exon_number}}=[$bChr,$bSt,$bEn];	
	}

	
}close(IN1);

close(IN1);
$NrAllSamples=scalar(keys %samplesH);
print STDERR "input features $c on ".scalar(keys %allowedChr)." reference sequences\n";


@aSamples=@{$$V{samples}};
print STDERR "samples in vcf: ".join(",",@aSamples)."\n";


print STDERR "open the bw files...\n";
############################################################
our %bwObj;

foreach $cSample (@aSamples){
	
	foreach $cT (("mq","q0","q20")){
				
# 		$cQBW=$jobfs."/".$cSample.".$cT.bw";
		$cQBW=$projectDir."/alignments/$cSample/bw/".$cSample.".$cT.bw";
		die " $cQBW ! exists" if (! -f $cQBW);
		$bwObj{$cSample}{$cT} = Bio::DB::Big->open($cQBW);
	}	
	print STDERR "# sample: $cSample, bw found\n";
}



#read in chr2len
open(IN1, "<$inRef.fai") || die "can not open $inRef";
while(<IN1>){ chomp; @_=split("\t",$_); $chr2len{$_[0]}=$_[1]; }close(IN1);


print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
############################################################
%meanCov;
foreach $cSample (@aSamples){
	$sumCov=0;
	for(1..22){
		$sumCov+=getCov($chrPf.$_,20000001,30000000,"q0",58,$cSample);
	}
	foreach $cChr ((keys %chr2len,"A")){ $meanCov{$cChr}{$cSample}=sprintf("%.2f",($sumCov)/22); }
	
	# determine X and Y coverage and round to 0, 0.5, 1, 1.5, 2, ... times average coverage
	$meanCov{$X_chr}{$cSample}=sprintf("%.0f",(getCov($X_chr,1,$chr2len{$X_chr},"q0",58,$cSample)/$meanCov{A}{$cSample})/0.5)*0.5*$meanCov{A}{$cSample};
	$meanCov{$Y_chr}{$cSample}=sprintf("%.5f", (getCov($Y_chr,6641419,10079253,"q0",0,$cSample)+getCov($Y_chr,13800704,23668908,"q0",0,$cSample))/2 ); # two mq 60 region one of each Y arm 
	print STDERR "Y:".$meanCov{$Y_chr}{$cSample}."\n";
	$meanCov{$Y_chr}{$cSample}=0 if $meanCov{$Y_chr}{$cSample} == (-1); # if -1 then because it is not covered	
	$meanCov{$Y_chr}{$cSample}=sprintf("%.0f",($meanCov{$Y_chr}{$cSample}/$meanCov{A}{$cSample})/0.5)*0.5*$meanCov{A}{$cSample};

	# MT coverage
	$meanCov{$M_chr}{$cSample}=sprintf("%.2f",(getCov($M_chr,1,$chr2len{$M_chr},"q0",58,$cSample)));	

	print STDERR "average coverage $cSample autosome:".$meanCov{A}{$cSample}.", X:".$meanCov{$X_chr}{$cSample}.", Y:".$meanCov{$Y_chr}{$cSample}.", MT:".$meanCov{$M_chr}{$cSample}."\n";
} 


print STDERR "annotating vcf...\n";
##################################################################
my $tabix = new Tabix(-data => "$inGenes.gff.gz",-index => "$inGenes.gff.gz.tbi");


$variantsAffectedByGene=0;
$allVariantsInVCF=0;
%var2bp2tr;
while ($V=$V->next_line){
  		
  	$allVariantsInVCF++;
  	undef %affectedGenesN;
	undef @affectedGenesD;
	undef %GTF;	
	
  	($cChr1,$bp1)=($$V{CHROM},$$V{POS});
  	
  	# next if no variant is at least PASS or high
  	$ph_count=0; foreach $cSample (@aSamples){  $ph_count++ if lc($$V{GTS}{FT}{$cSample}) eq "high" or lc($$V{GTS}{FT}{$cSample}) eq "pass";    }
  	next if $ph_count==0;
  	
  	$uniqid{$$V{"ID"}}++; die "id: ".$$V{"ID"}." is repeated!, needs to be unqiue otherwise variants get overwritten " if $uniqid{$$V{"ID"}}!=1;
#   	print STDERR ${$V->print_line} if $verbose;
	
	if ($$V{"INFO"}{"SVTYPE"} eq "BND"){
# 		next if exists($$V{"INFO"}{"SECONDARY"});
		if($$V{"ALT"} =~ /([\[\]])(.+):([0-9]+)([\[\]])N*$/){ # N[7:112692989[        [1:64854682[N
			($brkt1,$cChr2,$bp2)=($1,$2,$3);
			
			if(   $$V{"ALT"} =~ /^[ATGCN]/ and $brkt1 eq "[" ){ $bpOri1=1; $bpOri2=(-1);   } # N[7:112692989[      # ___| ----- |___    
			elsif($$V{"ALT"} !~ /^[ATGCN]/ and $brkt1 eq "]" ){ $bpOri1=(-1); $bpOri2=1;   } # ]7:112692989]N      # |___ ------ ___|
			elsif($$V{"ALT"} =~ /^[ATGCN]/ and $brkt1 eq "]" ){ $bpOri1=1; $bpOri2=1;      } # N]7:112692989]      # ___| -----  ___|    
			elsif($$V{"ALT"} !~ /^[ATGCN]/ and $brkt1 eq "[" ){ $bpOri1=(-1); $bpOri2=(-1);} # [7:112692989[N      # |___ -----  |___  
			
			if($cChr1 eq $cChr2){ 	
				$$V{"INFO"}{"LEN"}=abs($bp1-$bp2);
				$$V{"INFO"}{"END"}=$bp2;
			}
		}else{die $$V{"all"} }	
	}else{
		($cChr2,$bp2)=($$V{CHROM},$$V{"INFO"}{"END"});
		if($$V{"INFO"}{"SVTYPE"} eq "DEL"){ $bpOri1=1; $bpOri2=(-1); }
		elsif($$V{"INFO"}{"SVTYPE"} eq "INV"){ $bpOri1=1; $bpOri2=1; } # or $bpOri1=(-1); $bpOri2=(-1);
		elsif($$V{"INFO"}{"SVTYPE"} eq "DUP"){ $bpOri1=(-1); $bpOri2=1; }
		next;
	}
		

	######### if copy neutral only check at breakpoint, else annotate everything in between	
	if(  (exists($$V{"INFO"}{"CNV"}) and $$V{"INFO"}{"CNV"}==1) ){ # or  $$V{"INFO"}{"SVTYPE"} eq "INV" ){
		
		######### check region overlap
		checkReg($cChr1,$bp1,$cChr2,$bp2,$V);
		
		# and check breakpoint, too
		
	}else{
		
		######### check breakpoint
		checkBR($cChr1,$bp1,$bpOri1,$cChr2,$bp2,$bpOri2,$V);
		
	}
	
	
	
	# for each affected gene # %{$GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"}}}=%e;
	foreach $cGeneID (sort keys %GTF){
# 		print STDERR $cGeneID."\n";
		
		foreach $cTrID (  sort {  $tid2{exonic_len}{$b} <=> $tid2{exonic_len}{$a}       } keys %{$GTF{$cGeneID}}  ){ # sort by summed exon len
			
			
# 			print STDERR $cTrID."\n";
			undef %oTxt;
		
			###### gene name, gene ID, CANDG, biotype, transcript ID, affected feature type (breakpoint location CDS, exon, intron, intergenic) 
			$oTxt{gene_name}=$gid2{gene_name}{$cGeneID};
			$oTxt{gene_id}=$cGeneID;
			$oTxt{CANDG}=(exists( $geneList{$oTxt{gene_name}} ))? 1:0;
			$oTxt{biotype}=$gid2{gene_biotype}{$cGeneID};
			$oTxt{tr_id}=$cTrID;
			$oTxt{exonic_len}=$tid2{exonic_len}{$cTrID};
		
			foreach $cVarID (sort keys %{$GTF{$cGeneID}{$cTrID}} ){ # sort by summed exon len
				
				
				@bps=%{$GTF{$cGeneID}{$cTrID}{$cVarID}};
				die "error: don't expect more than two breakpoints per variant $cVarID" if @bps>2;
						
				
				@cOutSV=();
				%e=%{$GTF{$cGeneID}{$cTrID}{$cVarID}{$bps[0]}};
				# both_bp_in_gene:0, bp1_intronic:1, bp2_intronic:-1, cds_bases_affected:0, cds_bases_affected_pc:0, cds_bases_all:0,
				# cds_bases_unaffected:0, cds_len:, exonic_bases_affected:334, exonic_bases_affected_pc:88, exonic_bases_all:379,
				# exonic_bases_unaffected:45, exonic_len:379, exons_affected:1, exons_affected_txt:1-1, exons_h_affected:HASH(0x244cd660),
				# exons_h_unaffected:HASH(0x244cceb0), exons_unaffected:1,
				
				
				if($e{cds_bases_affected}>0){
					$oTxt{affected_feature}="CDS";
				}elsif($e{exonic_bases_affected}>0){
					$oTxt{affected_feature}="exon";
				}elsif($e{bp1_intronic} == 1 and $e{bp2_intronic} == 1){
					$oTxt{affected_feature}="intron";
				}else{
					$oTxt{affected_feature}="transcript";
				}
				
				map{ push @cOutSV, $oTxt{$_}  } @outCatGene;
				map{ push @cOutSV, $e{$_}  } @eCatGene;
				
				print STDERR "GTF_out $cGeneID $cTrID $cVarID $e{exons_affected_txt} ".join(" ", (map{ $e{$_}  } @eCatGene)  )."\n";
				
				
				# fusion partner, % affected CDS, inframe (yes,no), fused domain, known fusion cosmic,mittleman, 
				if($e{both_bp_in_gene}==0){ # $e{both_bp_in_gene}=0;
					
					
					$var2bp2tr{$$V{"ID"}}{"$bpChr:$bp"}{$cTrID}=$GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"}};
					
					$var2bp2tr{$$V{"ID"}}{"$bpChr:$bp"}
					
					$cTrID;
					

					
					%e=%{$GTF{$cGeneID}{$cTrID}{$cVarID}{$bps[1]};
					
					map{ push @cOutSV, $e{$_}  } ("gene_name",@eCatGene);
					
					
				}else{
					map{ push @cOutSV, "."  } ("gene_name",@eCatGene);
				}
				
				# get fusion partner in same line
				
				# get CNV
				
				
				###### some SV fields IA, IUA, IGV,GOTO
				map{ push @cOutSV, getVal($_)  } @annot_vcf_info;
				foreach $cSample (@aSamples){ map{ push @cOutSV, getVal($_)  } @annot_vcf_gts }
				
				
				print join("\t",(@cOutSV)  )."\n";
		
				
			}
			last; # only the longest transcript

		}
	
	}
	# 	last if $c2++>100;
	
	
	
}




sub checkBR{
		
		my ($cChr1,$bp1,$bpOri1,$cChr2,$bp2,$bpOri2,$V)=@_;
		@testCoords=([$cChr1,$bp1,$bpOri1,$cChr2,$bp2,$bpOri2],[$cChr2,$bp2,$bpOri2,$cChr1,$bp1,$bpOri1]);
		my $verbose=1;
		foreach $cA (@testCoords){
			
			my ($cChr,$cPos,$cOri,$oChr,$oPos,$oOri)=@$cA;
			
			if(!exists($allowedChr{$cChr}) ){return 0}
			print STDERR "--checkBR $cChr,$cPos,$oChr,$oPos\n" if $verbose;
		
			$res = $tabix->query($cChr,($cPos-5),($cPos+5)); 
			# 1. if only one breakpoint hits a gene, then it is affected by SV
			# 2. if both hit the gene but not the same intron, then the gene is also affected
			# in other words (since introns are not in gff), if overlap with any feature but the gene or transcript feature, then its affected
			if (defined $res->get){
			while(my $line = $tabix->read($res)){ #1       lincRNA exon    29554   30039   .       +       .       gene_id="ENSG00000243485";transcript_id="ENST00000473358";exon_number="1";gene_name="MIR1302-10";gene_source="ensembl_havana";gene_biotype="lincRNA";transcript_name="MIR1302-10-001";transcript_source="havana";exon_id="ENSE00001947070";
				
				chomp($line); ($bChr,$bSrc,$bFType,$bSt,$bEn,$bX1,$bOri,$bX3,$bName)=split("\t",$line);
				undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=substr($cTMP2[1],1,(length($cTMP2[1])-2)); }
				
				### only consider transcripts # all genes have a transcript!!!
				
				next if $bFType ne "transcript" or $bFType eq "Selenocysteine";
				die "!transcript_id $line\n" if !exists($h{transcript_id});
				die "!gene_id $line\n" if !exists($h{gene_id});
				die "!gene_name $line\n" if !exists($h{gene_name});
				
				# remaining: CDS,UTR,exon,start_codon,stop_codon
				$affectedGenesN{$h{gene_name}}++;
								
				# variant to breakpoint to affected transcript
				next if exists($var2bp2tr{$$V{"ID"}}{"$cChr:$cPos"});
				get_aff_exon_stat($cChr,$cPos,$cOri,\%h,$cChr1,$bp1,$bpOri1,$cChr2,$bp2,$bpOri2,$V);
				
				
			}}
		}
		
}


sub checkReg{
		
		my ($cChr1,$bp1,$cChr2,$bp2,$V)=@_;
		
		my $verbose=1;
		($bp1,$bp2)=($bp1<$bp2)? ($bp1,$bp2):($bp2,$bp1);
		print STDERR "--checkReg $cChr1,$bp1,$cChr2,$bp2\n" if $verbose;
		
		my ($qChr,$qPos1,$qPos2)=($cChr1,($bp1-5),($bp2+5));
		
		if(!exists($allowedChr{$qChr}) ){return 0}
		
		$res = $tabix->query($qChr,$qPos1,$qPos2); 
	
		# 1. if any of these features "CDS,UTR,exon,start_codon,stop_codon" are hit, the gene is affected
		while(my $line = $tabix->read($res)){ #1       lincRNA exon    29554   30039   .       +       .       gene_id="ENSG00000243485";transcript_id="ENST00000473358";exon_number="1";gene_name="MIR1302-10";gene_source="ensembl_havana";gene_biotype="lincRNA";transcript_name="MIR1302-10-001";transcript_source="havana";exon_id="ENSE00001947070";
			
			chomp($line); ($bChr,$bSrc,$bFType,$bSt,$bEn,$bX1,$bOri,$bX3,$bName)=split("\t",$line);
			undef %h; foreach $cTMP (split(";",$bName)){ @cTMP2=split("=",$cTMP); $h{$cTMP2[0]}=substr($cTMP2[1],1,(length($cTMP2[1])-2)); }
			
			next if $bFType ne "transcript" or $bFType eq "Selenocysteine";
			die "!transcript_id $line\n" if !exists($h{transcript_id});
			die "!gene_id $line\n" if !exists($h{gene_id});
			die "!gene_name $line\n" if !exists($h{gene_name});
					
			# remaining: CDS,UTR,exon,gene,start_codon,stop_codon

			$affectedGenesN{$h{gene_name}}++;

			# variant to breakpoint to affected transcript
			next if exists($var2bp2tr{$$V{"ID"}}{"$cChr1:$bp1:$bp2"});
			get_aff_exon_stat($cChr1,$bp1,0,\%h,$cChr1,$bp1,0,$cChr2,$bp2,0,$V); # breakpoint orientation not needed if on same chr

			
		}
}


	
sub get_aff_exon_stat {
	
	my($bpChr,$bp,$bpOri,$h,$bpChr1,$bp1,$bpOri1,$bpChr2,$bp2,$bpOri2,$V)=@_;
	$cTrID=$$h{transcript_id};
	
	my $verbose=0;
	undef %e;
	
	if($bpChr1 eq $bpChr2 and $bSt>= $bp1 and $bp1<=$bEn and $bSt>= $bp2 and $bp2<= $bEn ){ # variant st and end is within gene
		$e{both_bp_in_gene}=1;
	}else{
		$e{both_bp_in_gene}=0;
	}
				
	die "VCF coordinates are not in the correct order: $bpChr1 eq $bpChr2 and $bp1 > $bp2" if $bpChr1 eq $bpChr2 and $bp1 > $bp2;
	$cTrH=$tid2{exons}{$cTrID};
	$trOri=$tid2{orentation}{$cTrID};	
	$transcript_st=$tid2{transcript_st}{$cTrID};
	$transcript_en=$tid2{transcript_en}{$cTrID};
	$trChr=$tid2{'chr'}{$cTrID};
	
	$e{gene_name}=$$h{gene_name};
	$e{transcript_id}=$$h{transcript_id};
			
	$e{exons_affected}=0;
	$e{exons_unaffected}=0;
	$e{exonic_bases_affected}=0;
	$e{exonic_bases_unaffected}=0;
	$e{exonic_bases_all}=0;
	$e{exonic_len}=$tid2{exonic_len}{$cTrID};
	
	$e{cds_bases_affected}=0;
	$e{cds_bases_unaffected}=0;
	$e{cds_bases_all}=0;
	$e{cds_len}=$tid2{cds_len}{$cTrID};
	
	$e{bp1_intronic}=(-1);
	$e{bp2_intronic}=(-1);
	# if st and end on diff chroms - one end	
	
	$e{exons_h_affected}={};
	$e{exons_h_unaffected}={};
	
	if($bpChr1 ne $bpChr2 ){ #if variant st and end are on different chromosomes
		$affOrUnaff=($bpOri*$trOri>0)? "unaffected":"affected";
	}elsif($bpChr1 eq $bpChr2 and $bp1<=$transcript_st and $transcript_st<=$bp2){ # variant overlaps gene start -> affected
		$affOrUnaff="affected";
	}elsif($bpChr1 eq $bpChr2 and ($transcript_st*$trOri)<($bp1*$trOri) and ($transcript_st*$trOri)<($bp2*$trOri) ){ # variant downstream of start, either within gene or partially outisde gene
		$affOrUnaff="unaffected";
	}else{
		die "variant versus gene location not accounted for chr: $bpChr1 ne $bpChr2 bp1:$bp1, bp2:$bp2, trOri:$trOri, transcript_st: $transcript_st \n";
	}	
	$iffOrUnaff=$affOrUnaff;
	
	$switched1=0;$switched2=0;
	print STDERR "bp:$bp, bpOri:$bpOri, cTrID:$cTrID, trOri:$trOri, affOrUnaff: $affOrUnaff ".($bpOri*$trOri>0)."\n" if $verbose;
	
	foreach $cExon ( sort {$a <=> $b} keys  %{$cTrH}  ){ # for each exon
		($eChr,$eSt,$eEn)=@{$$cTrH{$cExon}};
		($eSt,$eEn)=($trOri == (-1))? ($eEn*(-1),$eSt*(-1)):($eSt,$eEn);
		$res_test_br_in_exon=0;
		for($i=$eSt; $i<=$eEn; $i++ ){ # for each exon position
			$e{exonic_bases_all}++;
# 			print STDERR "i:$i >= ".($bp1*$trOri)."  ".($bp2*$trOri)." ";
			
			if( $i>=($bp1*$trOri) and $bpChr1 eq $trChr and $switched1==0){
				$affOrUnaff=($affOrUnaff eq "affected")? "unaffected":"affected";
				if( $i!=($bp1*$trOri) ){ $e{bp1_intronic}=1; }else{$res_test_br_in_exon=1}
				$switched1++;
				print STDERR "switch1\n" if $verbose;
			}
			if( $i>=($bp2*$trOri) and $bpChr2 eq $trChr and $switched2==0){
				$affOrUnaff=($affOrUnaff eq "affected")? "unaffected":"affected";
				if( $i!=($bp2*$trOri) ){ $e{bp1_intronic}=1; }else{$res_test_br_in_exon=1}
				$switched2++;
				print STDERR "switch2\n" if $verbose;
			}
			$e{"exonic_bases_".$affOrUnaff}++;
		}
		
		$lffOrUnaff=(  $res_test_br_in_exon )? "affected":$affOrUnaff;
		$e{"exons_h_".$lffOrUnaff}{$cExon}++;
		$e{"exons_".$lffOrUnaff}++;
		print STDERR "cExon:$cExon, eChr:$eChr, eSt:$eSt, eEn:$eEn,  affOrUnaff: $lffOrUnaff, res_test_br_in_exon:$res_test_br_in_exon \n" if $verbose;
	}
	
	$affOrUnaff=$iffOrUnaff;
	$switched1=0;$switched2=0;
	foreach $cExon ( sort {$a <=> $b} keys  %{$cTrH}  ){
		next if !exists($tid2{cds}{$cTrID}{$cExon});
		($cdsChr,$cdsSt,$cdsEn)=@{$tid2{cds}{$cTrID}{$cExon}};
		
		($eSt,$eEn)=($trOri == (-1))? ($cdsEn*(-1),$cdsSt*(-1)):($cdsSt,$cdsEn);
		for(my $i=$eSt; $i<=$eEn; $i++ ){
			$e{cds_bases_all}++;
			if( $i>=($bp1*$trOri) and $bpChr1 eq $trChr and $switched1==0){
				$affOrUnaff=($affOrUnaff eq "affected")? "unaffected":"affected";
				$switched1++;
				print STDERR "switch1\n" if $verbose;
			}
			if( $i>=($bp2*$trOri) and $bpChr2 eq $trChr and $switched2==0){
				$affOrUnaff=($affOrUnaff eq "affected")? "unaffected":"affected";
				$switched2++;
				print STDERR "switch2\n" if $verbose;
			}
			$e{"cds_bases_".$affOrUnaff}++;
		}
	}
	
	$eLen=$tid2{exonic_len}{$cTrID};
	die "exon lencht discrepancy $eLen != ".$tid2{exonic_len}{$cTrID} if $eLen != $tid2{exonic_len}{$cTrID};
	$e{exonic_bases_affected_pc}=round($e{exonic_bases_affected}/$e{exonic_bases_all}*100,0);
	$e{cds_bases_affected_pc}=($e{cds_bases_all}>0)? round($e{cds_bases_affected}/$e{cds_bases_all}*100,0):0;
	
	### get affected st and end exon
	@eTmp= sort { $a <=> $b } keys %{$e{exons_h_affected}};
	if(scalar(@eTmp)==0){ $e{exons_affected_txt}="NA"; }
	else{
		$ss=$eTmp[0];
		for($j=1;$j<=$#eTmp;$j++){
			if($eTmp[$j]>$eTmp[$j-1]+1){  $ss.="-$eTmp[$j-1],$eTmp[$j]" }
		}
		$ss.="-".$eTmp[$#eTmp];
	}
	$e{exons_affected_txt}=$ss;
	print STDERR "   ".scalar( keys %{$e{exons_h_affected}} )." exons_h_affected  $ss \n" if $verbose;
	
	die "error: second breakpoint for variant $$V{ID} \n" if exists($GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"});
	%{$GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"}}}=%e;
	
	$var2bp2tr{$$V{"ID"}}{"$bpChr:$bp"}{$cTrID}=$GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"}};
	
	if($verbose){ map {print STDERR "$_:".$e{$_}.", " } sort keys %e; }
	
# 	print STDERR "GTF:".$GTF{$$h{gene_id}}{$cTrID}{$$V{"ID"}}{exons_affected_txt}."  \n";# if $verbose;
# 	$a=<STDIN>;
	
	print STDERR "details hit BR:  ".scalar(keys %GTF)." $bFType ".$$h{gene_name}." ".$cTrID." $bChr:$bSt-$bEn geneOri:$bOri\n" if $verbose;
	
	
	return;
	
}


print STDERR "of all variants in VCF $allVariantsInVCF, $variantsAffectedByGene (".sprintf("%.1f",($variantsAffectedByGene/$allVariantsInVCF*100))."%) hit ENSEMBL genes\n";



sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}



sub getCov{
	my ($qChr,$qSt,$qEn,$covQ,$qMQ,$cSample)=@_;
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	
	$qSt=0 if $qSt<0;
	$qEn=$chr2len{$qChr} if $qEn > $chr2len{$qChr};
	return (-1) if $qSt>=$qEn;
	
	my $qLen=($qEn-$qSt+1);
	my ($sumCov,$fBins)=(0,0);
    my $qBins=100;

    $qBins=$qLen if ($qLen/$qBins)<1;

	my $aCov_stat = $bwObj{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	my $aMQ_stat = $bwObj{mq}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean') if $qMQ>0;	
	
	for ($i=0; $i<=$#$aCov_stat;$i++){		
		
      $cMQ=$qMQ>0 ? ${$aMQ_stat}[$i]:0;
      $cCov=${$aCov_stat}[$i];

	  next if $qMQ>0 and $cMQ<$qMQ;
	  next if $cCov<0;	  
	  $sumCov+=$cCov; $fBins++;
	  		
	}

	
	if($fBins==0 or $fBins/$qBins<0.33333){ return (-1) }
	else{ return sprintf("%.2f",$sumCov/$fBins) }
}



sub getVal {

	($cCol)=@_;

	if($cCol eq "PHEN"){
		undef @tmpTxt;
		$cGeneList=$$V{INFO}{GENES};
		$cGeneList=~ s/"//g;
		foreach $cGene (split(",",$cGeneList)){
			if ( exists($ph_gene2phen{$cGene}) ){
				push @tmpTxt, "$cGene,".$ph_gene2phen{$cGene};
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
			($cChr2,$bp2)=($1,$2);
			}
		}else{
			($cChr2,$bp2)=($$V{CHROM},$$V{"INFO"}{"END"});
		}
		
		
		$flan2=1000;
		
		if( $$V{CHROM} eq $cChr2 and $$V{INFO}{CNV}==1 ){ $flan=int(($$V{INFO}{"END"}-$$V{POS})*1.5);  $min1=($$V{POS}-$flan); $min1=1 if $min1<=0;  $cReg=$$V{CHROM}.":".$min1."-".($$V{INFO}{"END"}+$flan)   }
		else{  $min1=($$V{POS}-$flan2); $min1=1 if $min1<=0; $min2=($bp2-$flan2); $min2=0 if $min2<0;  $cReg=$$V{CHROM}.":".$min1."-".($$V{POS}+$flan2)."%20".$cChr2.":".$min2."-".($bp2+$flan2) }
		
		return "=HYPERLINK(\"http://localhost:60151/goto?locus=$cReg\",\"GOTO\")";
		
	}elsif($cCol eq "HPO"){
		
		if(exists($$V{INFO}{$cCol}) and $$V{INFO}{$cCol} !~ /^[|]+$/){ return $$V{INFO}{$cCol} }
		else{ return "" }
	}
	
	
	
	return $cSample if $cCol eq "SAMPLE";
	
	return $sa2info{$cSample}{$cCol} if exists($sa2info{$cSample}{$cCol});
	return $$V{GTS}{$cCol}{$cSample} if exists($$V{GTS}{$cCol}{$cSample});
	return $$V{INFO}{$cCol} if exists($$V{INFO}{$cCol});	
	return $$V{$cCol} if exists($$V{$cCol});	
	return ".";
}



sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}

