






use Getopt::Long;



GetOptions ("gene=s" => \$geneList,
			"hpo=s"  => \$hpoList,
		    "cyto=s"  => \$cytoList,
		    "reg=s"  => \$regList,
		    "v"  => \$verbose
)  or die("Error in command line arguments\n");




# read ENSEMBL gene IDS







if($geneList){
	open(IN1, "</g/data2/gd7/research/andmin/resources/genes/Homo_sapiens.GRCh37.75.ids") || die " nicht gefunden";
	while(<IN1>){ chomp; 
	$validGeneIDs{$_}++;
	}close(IN1);

	if($geneList =~ /\,/ or   $geneList !~ /\//  ){
		map { $inclGenes{$_}++; } split(",",$geneList);
		
	}else{
		open(IN1, "<$geneList") || die "$geneList not found"; 
		while(<IN1>){ chomp; $_=~ s/^ +| +$//g; $inclGenes{$_}++;  }close(IN1); 
	}
	print STDERR scalar(keys %inclGenes)." genes\n";
	map { print STDERR "gene name $_ ! exists \n" if !exists($validGeneIDs{$_}) } keys %inclGenes;
}
print STDERR "done checking genes\n";

if($hpoList){
	if($hpoList =~ /\,/ or   $hpoList !~ /\//  ){
		map { $inclHPO{$_}++ } split(":",$hpoList);
	}else{
		open(IN1, "<$hpoList") || die "$hpoList not found"; 
		while(<IN1>){ chomp; $inclHPO{$_}++; }close(IN1); 
	}
	print STDERR scalar(keys %inclHPO)." HPOs\n";
}

if($regList){
	%inclReg;$inclRegC=0;$inclRegS="";
	foreach $cR (split(",",$regList)){
		$cR=~ s/,//g;
		@tmp1=split(":",$cR);
		@tmp2=split("-",$tmp1[1]);
		$inclRegS.="$tmp1[0]:$tmp2[0]-$tmp2[1], ";
		push @{$inclReg{$tmp1[0]}}, [$tmp2[0],$tmp2[1]]; $inclRegC++;
	}
	print STDERR $inclRegC." regions: $inclRegS\n";
}

if($cytoList){ #11q12.2-11q12.3,16q22.1

	open(IN1, "</g/data2/gd7/resources/UCSC/GRCh37/cytoBand.txt") || die " nicht gefunden";
	while(<IN1>){ chomp; @_=split("\t",$_); $C++; $cyto2C{$_[0].$_[3]}=$C; $C2cyto{$C}=$_[0].$_[3]; $cyto2Loc{$_[0].$_[3]}=[$_[0],$_[1],$_[2]];
	}close(IN1);
	
	%inclReg;$inclRegC=0;
	if($cytoList =~ /\,/ or   $cytoList !~ /\//  ){
		map { readCyto($_) } split(",",$cytoList);
	}else{
		open(IN1, "<$cytoList") || die "$cytoList not found"; 
		while(<IN1>){ chomp; readCyto($_) }close(IN1); 
	}
	print STDERR $inclRegC." cytband regions\n";
}

sub readCyto {
	
	$sCyt=shift;
	undef @Cyts;
	
	if ($sCyt =~ /-/){		
		($cytSt,$cytEn)=split(/-/,$sCyt);
		die "cytoband $cytSt does not exists!" if !exists($cyto2C{$cytSt});
		die "cytoband $cytEn does not exists!" if !exists($cyto2C{$cytEn});
		for $C ($cyto2C{$cytSt}..$cyto2C{$cytEn}){  push @Cyts,  $C2cyto{$C};}	
	}else{
		push @Cyts, $sCyt;
	}	
	foreach $cCyt (@Cyts){
		die "cytoband $cCyt does not exists!" if !exists($cyto2Loc{$cCyt});
		($chr,$st,$en)=@{$cyto2Loc{$cCyt}};
		push @{$inclReg{$chr}}, [$st,$en]; $inclRegC++;
		print STDERR "use cytoband $cCyt $chr:$st-$en\n";	
	}
}



while(<STDIN>){

	if(/^#/){ print; next } 

	$printOK=0;

	if ($geneList and /GENES="([^;]+)";/){
		@tmp=split(",",$1);
		foreach $cGene (@tmp){ $printOK++ if exists($inclGenes{$cGene});   } 
	}

	if ($hpoList and /HPO="([^;]+)";/){
		@tmp=split(/[,|]/,$1);
		foreach $cHPO (@tmp){ $printOK++ if exists($inclHPO{$cHPO});   } 
	}
	
	
	if(%inclReg){
		undef(@testPos);
		@_=split("\t",$_);
		push @testPos, $_[1];
		if(/END=([0-9]+);/){ push @testPos, $1 }
		
		foreach $cPos (@testPos){
			foreach $cA (@{$inclReg{$_[0]}}){
				($st,$en)=@$cA;
				if ($st<=$cPos and $cPos<=$en){ $printOK++ }
			}
		}
	}
	
	
	print if $printOK;

}





sub overlap{ my($R1,$R2,$Q1,$Q2)=@_; ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;}















