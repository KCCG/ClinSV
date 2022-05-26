
# Author: Andre Minoche, a.minoche@garvan.org.au


%allowedChr; for(1..22){ $allowedChr{$_}++ } 
$allowedChr{"X"}++; $allowedChr{"Y"}++;
 
%allowedChr; for(1..22){ $allowedChr{"chr".$_}++ } 
$allowedChr{"chrX"}++; $allowedChr{"chrY"}++; 


$maxLinkSize=500000;
$readLength=150;
$fl=1000;

$cMode=shift(@ARGV);
$outStem=shift(@ARGV);
$cSample=shift(@ARGV);
$ref2=shift(@ARGV);

@inFileS=@ARGV;
$verbose=0;

%validDPIDs;

@greyCols=("#362e2f","#373c40","#394543","#464646","#4e5159","#564c56","#6a5c5d","#727272","#7d7b78","#988f89","#918582","#7e7a6d","#7d7c67","#7f7479","#677479","#969696","#a7a0a1","#9e9c90","#9ea3a1","#bebebe","#bbbbb3","#bdb9b5","#bcbfbd","#cacacc");

print STDERR "input bams: ".join(",",@inFileS).";\n";

# 		inversion:= lila 169:104:181, a968b5
# 		dup:= blue = 0375b6 = 3,117,182
# 		del:= red = c34343 =195,67,67
# 		BND:= grey = 5d5553 = 93,85,83

#		orphan: yellow = f8a748 = 248,167,72 # orphans are filtered out

%PE_type2col=("DEL","#c34343","DUP","#0375b6","INV","#a968b5","BND","#5d5553");

%SR_type2col;%SR_type2name;
map { $SR_type2col{$_}="#c34343"; $SR_type2name{$_}="DEL"; } ("1pU,2pD","2mU,1mD"); # del type
map { $SR_type2col{$_}="#0375b6"; $SR_type2name{$_}="DUP"; } ("1mU,2mD","2pU,1pD"); # dupe type
map { $SR_type2col{$_}="#a968b5"; $SR_type2name{$_}="INV"; } ("2mU,1pD","2pU,1mD","1pU,2mD","1mU,2pD"); # inv type
map { $SR_type2col{$_}="#5d5553"; $SR_type2name{$_}="BND"; } ("1p,2p","1m,2m","1m,2p","1p,2m"); # different chromosomes
# map { $SR_type2col{$_}="#5d5553"; $SR_type2name{$_}="ORPH"; } ("2m","2p","1p","1m"); # orphan type


sub type2colFun{

	my ($PEorSR,$abrType,$xchr)=@_;
	$chr_= $xchr =~ /^chr/ ? substr($xchr,3):$xchr;
	$cGreyCol=$greyCols[23];
	if($chr_ =~ /^[0-9]+$/){ $cGreyCol=$greyCols[(($chr_*1)-1)]; }
	elsif($xchr =~ /^X$/){ $cGreyCol=$greyCols[22]; }

	if($PEorSR eq "SR"){
		if( $SR_type2name{$abrType} eq "BND"){ return $cGreyCol}
		else{ return $SR_type2col{$abrType} }
	}elsif($PEorSR eq "PE"){
		if( $abrType eq "BND"){ return $cGreyCol}
		else{ return $PE_type2col{$abrType} }
	}

}

@aCluster;

open(IGV, ">$outStem/bed/$cSample.$cMode.bed") || die " can not write to $outStem/bed/$cSample.$cMode.bed";
open(BRKP, ">$outStem/bed/$cSample.$cMode.brkp") || die " can not write to $outStem/bed/$cSample.$cMode.brkp";
open(SAM, " | samtools view -b - > $outStem/$cSample.$cMode.f.bam") || die " can not write to $outStem/$cSample.$cMode.f.bam";
print IGV "#gffTags\n";

####### print header needed for multiple discordant bam files
for ($q=0; $q<=$#inFileS;$q++){

	$inFile=$inFileS[$q];
	print STDERR "reading header $inFile\n";
	
	open(IN1, "samtools view -H $inFile | ") || die "  $inFile nicht gefunden";
	while(<IN1>){
		if ($q==0){print SAM $_; next}
		elsif (/^[@]RG/){print SAM $_; next}
	}
	close(IN1);
}


%hPosID;

for ($q=0; $q<=$#inFileS;$q++){

	$inFile=$inFileS[$q];
	
	exit 1 if (! -e $inFile);
	
	print STDERR "reading $inFile\n";
	
	open(IN1, "samtools view $inFile | ") || die "  $inFile nicht gefunden";
	while(<IN1>){
	
# 		if (/^@/){next}
	#     last if $c++<100;
		chomp; @s=split("\t",$_);
		$stat{all_alignments}++;
		if ($s[1] & 4){ $stat{read_unmapped}++; print "read_unmapped $_\n" if $verbose; next;}
		if ($s[1] & 8){ $stat{mate_unmapped}++; print "mate_unmapped $_\n" if $verbose; next;}
		$subID=0;
		$printR1=1;
		($m1Chr,$m1St,$m1Cigar,$cID,$m1MQ)=($s[2],$s[3],$s[5],readID($s[0]),$s[4]); 
		($m2Chr,$m2St)=($s[6],$s[7]);
		$m2Chr=$m1Chr if ($m2Chr eq "=");
		
		$m1En=($m1St-1+findEnd($m1Cigar));
		$m1Len=abs($m1En-$m1St)+1;
		$m1Ori=($s[1] & 16)? "-":"+";
	# 	last if $m1Chr eq "2";
		
		if (abs($m1St-$m2St)<=5){ $stat{read1_and_read2_same_start}++; next;}
		
		if ($s[1] & 1024){ $stat{pcr_duplicate}++; print "pcr_duplicate $_\n" if $verbose; next;}
	
		if (clipPartCount($m1Cigar)>=2){$stat{more_than_two_clipped_parts}++; print "more_than_two_clipped_parts $_\n" if $verbose; next;} # these reads are probably not present in the reference. Middle parts of 36 bases are too unspecific and often match to several locations, their partner should be considered orphan
		# in STR rerions some split reads are clipped more than once, this can be exploited further.	
		#next if $m1MQ<20;
	
	# 	print  "p1 $m1Chr:$m1St\-$m1En, $cID, $m1MQ, $m1Ori $s[5] F:$s[1] $isPrimAlign ".clipPartCount($m1Cigar)."\n";
 
		############ output split reads ############  Flag SA:Z:5
		if(substr($cMode,0,1) eq "s" ){
			#print " $_\n";
			$m1MLoc=getMatchLocationOnRead($m1Ori,$m1Cigar); # returns match length. If positive at the begging of read if negative at the end of the read
			if ($m1MLoc==0){
				$stat{no_clipped_part}++; 
				print "no_clipped_part $_\n" if $verbose;
				next;
			}
					
			# get current SA:Z string, next if none is present (no split read)
			$cSAZ=""; for $cF (11..$#s){ if (substr($s[$cF],0,4) eq "SA:Z"){$cSAZ=$s[$cF];last} }
			if ($cSAZ eq ""){ 
				#my $m1Type=($m1MLoc>0)? "2":"1"; $m1Type.=($m1Ori eq "+")? "p":"m";
				#$gffTags="Name=".$SR_type2name{$m1Type}."_$m1Type;readID=$cID;cigar=$m1Ori:$m1Cigar;";
				#print IGV join("\t",($m1Chr,($m1St-1),$m1En,$gffTags,$m1MQ,$m1Ori,($m1St-1),$m1En,type2colFun("SR",$m1Type,$m2Chr),2,"$m1Len","0"))."\n";
				$stat{no_SA_string}++; print "no_SA_string $_\n" if $verbose;
				next;
			} # print orphan split read
				
			# only take in to account split reads with two aligning portions
			@splitMaps=split(";",$cSAZ);
			if (scalar(@splitMaps)>1){$stat{more_than_two_SA_alignments}++; print "more_than_two_SA_alignments $_\n" if $verbose; next;}
		
			# extract the coordinates from $cSAZ
			# SA:Z:HLA-B*55:24,3004,-,20S129M,0,3;
			
			($sChr,$sSt,$sOri,$sCigar,$sMQ)=split(",",substr($splitMaps[0],5)); 
			$sEn=$sSt-1+findEnd($sCigar); 
			$sLen=abs($sEn-$sSt)+1; 
			if (clipPartCount($sCigar)>=2){
				$stat{more_than_two_clipped_parts}++; 
				print "more_than_two_clipped_parts $_\n" if $verbose;
				next;
			}
			
			#next if $sMQ<20;
			
			# next if same chr,start,cigar for part1 and part2 already exists. 
			$cPosID="$m1Chr,$m1St,$m1En,$sChr,$sSt,$sEn";$cPosID2="$sChr,$sSt,$sEn,$m1Chr,$m1St,$m1En";	
			if (  (exists($hPosID{$cPosID}) and $hPosID{$cPosID} ne $s[0]) or (exists($hPosID{$cPosID2}) and $hPosID{$cPosID2} ne $s[0]) ){
				$stat{same_split_read_coordinates_already_exists}++; print "same_split_read_coordinates_already_exists xxxxxxx $_\n" if $verbose; 
				next;
			}else{ $hPosID{$cPosID}=$s[0]; }
			
			
			# which is first part of the read which is the second part of the read
			# orientate within read using the soft and hard clipped parts, because matching parts can be interrupted by indels
			$sMLoc=getMatchLocationOnRead($sOri,$sCigar);
			if ($sMLoc==0){
				$stat{no_clipped_part}++; 
				print "no_clipped_part $_\n" if $verbose;
				next;
			}
			
			if (($m1MLoc<0 and $sMLoc<0) or ($m1MLoc>0 and $sMLoc>0)){$stat{parts_on_read_overlap}++; print "parts_on_read_overlap $_\n" if $verbose;  next;} # different parts of the read have to match, else it is no split read
			
			
			# the different parts of the split read must be either close to read1 or read2, else this would be odd
			if(   ($m1Chr ne $sChr or abs($sSt-$m1St)>2000 )  and  ($m2Chr ne $sChr or abs($sSt-$m2St)>2000 ) 	and  ($m1Chr ne $m2Chr or abs($m1St-$m2St)>2000 )){
				#print "# split read far from read1 and read2 $m1Chr ne $sChr $sSt-$m1St && $m2Chr ne $sChr $sSt-$m2St && $m1Chr ne $m2Chr or abs($m1St-$m2St) $_\n";
				$stat{parts_mapping_at_differnet_locations}++; print "parts_mapping_at_differnet_locations $_\n" if $verbose; 
				next;
			}
			
			$stat{passing_alignments}++;
			print SAM $_."\n";
			
			# The secondary alignment corresponds to the second match part for split reads. Do not exclude these of the bam, since this is what lumpy actually uses, not the SA tags, as I do. 
			if ($s[1] & 256){ next;}
			
			
			$m1Dist=($sSt-$m1St); # distance between mapping 1 start and mapping 2 start, if positive this mapping is upstream
			$sDist=($m1St-$sSt);
		
			($m1Type,$sType,$mType)=getSRType($m1Chr,$m1MLoc,$m1Ori,$m1Dist,$sChr,$sMLoc,$sOri,$sDist);
			print "$m1Chr:$m1St-$m1En:$m1Ori:$m1Cigar:L$m1MLoc:D$m1Dist:$m1Type | $sChr:$sSt-$sEn:$sOri:$sCigar:L$sMLoc:D$sDist:$sType \n" if !exists($SR_type2name{$mType});
			
			
			#### what are the breakpoint coordinates, make the breakpoint side appear larger in IGV   ---#   #---- or #---   ---#, # thickEnd only accepts one coord, so only possible for DEL and non connected features
			if($m1Type =~ /^1p/ or $m1Type =~ /^2m/){   $m1Brkp=$m1En; $m1NBrkp=$m1St;  $m1ThickStEn=($m1En-5-1)."\t".$m1En;  }
			elsif($m1Type =~ /^1m/ or $m1Type =~ /^2p/){   $m1Brkp=$m1St; $m1NBrkp=$m1En; $m1ThickStEn=($m1St-1)."\t".($m1St+5); }
			
			if($sType  =~ /^1p/ or $sType  =~ /^2m/ ){   $sBrkp=$sEn; $sNBrkp=$sSt;  $sThickStEn=($sEn-5-1)."\t".$sEn;  }
			elsif($sType  =~ /^1m/ or $sType =~ /^2p/){   $sBrkp=$sSt; $sNBrkp=$sEn;  $sThickStEn=($sSt-1)."\t".($sSt+5); }
			
			# sort breakpoint location now and at retrieval the same way to prevent the need to search for each breakpoint in case on different chromosomes
			@rp = sort { $$a[0] cmp $$b[0] || $$a[1] <=> $$b[1];  } ([$m1Chr,$m1Brkp,$m1Type,$m1NBrkp,$m1MQ],[$sChr,$sBrkp,$sType,$sNBrkp,$sMQ]);
			print BRKP join("\t",($rp[0][0],$rp[0][1],"type=".join(",",( sort ($rp[0][2],$rp[1][2]) )).";type1=$rp[0][2];readID=$cID;NBP1=$rp[0][3];MQ1=$rp[0][4];chr2=$rp[1][0];BP2=$rp[1][1];type2=$rp[1][2];NBP2=$rp[1][3];MQ2=$rp[1][4]" ) )."\n";		
			
					
			
			###### output as bed for IGV		
			if ($m1Chr eq $sChr and abs($m1Dist)<$maxLinkSize ){
			
# 				next if !exists($allowedChr{$m1Chr});
			
				$gffTags="Name2=".$SR_type2name{$mType}."_$mType;readID=$cID;MQ1=$m1MQ;MQ2=$sMQ;distance=".abs($m1Dist).";cigar1=$m1Ori:$m1Cigar;cigar2=$sOri:$sCigar;";	
				#$gffTags=$mType\_".abs($m1Dist)."_".$cID";
				if($m1St<=$sSt){
					print IGV join("\t",($m1Chr,($m1St-1),$sEn,$gffTags,$m1MQ,$m1Ori,($m1St-1),$sEn,type2colFun("SR",$mType,$sChr),2,"$m1Len,$sLen","0,$m1Dist"))."\n";
				}else{
					print IGV join("\t",($m1Chr,($sSt-1),$m1En,$gffTags,$sMQ,$sOri,($sSt-1),$m1En,type2colFun("SR",$mType,$sChr),2,"$sLen,$m1Len","0,$sDist"))."\n";		
				}
			}else{
				$gffTags="Name2=".$SR_type2name{$mType}."_$m1Type;readID=$cID;MQ=$m1MQ;loc2=$sChr:$sSt;cigar=$m1Ori:$m1Cigar;joinType=$mType;";
				print IGV join("\t",($m1Chr,($m1St-1),$m1En,$gffTags,$m1MQ,$m1Ori,$m1ThickStEn,type2colFun("SR",$mType,$sChr),2,"$m1Len","0"))."\n";# if exists($allowedChr{$m1Chr});
			
				$gffTags="Name2=".$SR_type2name{$mType}."_$sType;readID=$cID;MQ=$sMQ;loc2=$m1Chr:$m1St;cigar=$sOri:$sCigar;joinType=$mType;";
				print IGV join("\t",($sChr,($sSt-1),$sEn,$gffTags,$sMQ,$sOri,$sThickStEn,type2colFun("SR",$mType,$m1Chr),2,"$sLen","0"))."\n";# if exists($allowedChr{$sChr});
			}
			
		}
	
		############ 
		############ discordant reads
		if(substr($cMode,0,1) eq "d" ){
		
		
			$isPrimAlign=($s[1] & 256)? 0:1;
			if (!$isPrimAlign){ $stat{not_primary_alignment}++; next;}
	# 		next if $m1MQ<20;
	# 		next if !exists($allowedChr{$m1Chr});
			
			if ($s[1] & 64){ $stat{not_first_read_in_pair}++; next;}
		
			$readReverseStrand=($s[1] & 16)? 1:0;
			$mateReverseStrand=($s[1] & 32)? 1:0;
					
			
	# 		next if !exists($allowedChr{$m2Chr});
		
			# get the end pos of mate, MC:Z: and MQ of mate
			$MCp=0;
			for $cF (11..$#s){  
				if (substr($s[$cF],0,4) eq "MC:Z"){  ($x1,$x2,$m2Cigar)=split(":",$s[$cF]);  $m2En=$m2St-1+findEnd($m2Cigar); $m2Len=abs($m2En-$m2St)+1; $MCp++}
				if (substr($s[$cF],0,4) eq "MQ:i"){  ($x1,$x2,$m2MQ)=split(":",$s[$cF]); }
			}
		
			if ($MCp==0){
				$m2Cigar=length($s[9])."M";
				$m2En=$m2St-1+findEnd($m2Cigar); $m2Len=abs($m2En-$m2St)+1;
				$m2MQ="-1";
			}else{
	# 			next if $m2MQ<20;
			}
		
	# 		print  "no MC flag present run samblaster with --addMateTags !?\n" if $MCp==0;
		
			#print  STDERR "p1 $m1Chr:$m1St\-$m1En, $cID, $m1MQ, $m1Ori $s[5] F:$s[1] $isPrimAlign\n";
			#print  STDERR "p2 $m2Chr:$m2St\-$m2En, $cID, $m2MQ, $m2Ori \n";
		
			@AMM=sort {$a<=>$b} ($m1St,$m1En,$m2St,$m2En); 
			($cMin,$cMax)=($AMM[0],$AMM[$#AMM]);
			$cGlobLen=($cMax-$cMin+1);
			$m1Ori=($readReverseStrand==1)? "-":"+";
			$m2Ori=($mateReverseStrand==1)? "-":"+";
			
			$m1Type=($readReverseStrand==1)? "m":"p";
			$m2Type=($mateReverseStrand==1)? "m":"p";
			
			if($s[3]<=$s[7]){ 
				($UpRev,$DoRev)=($readReverseStrand,$mateReverseStrand);
				@ALen=($m1Len,$m2Len);
				@ASt=(($m1St-$cMin),($m2St-$cMin));
				
				if($m1Chr eq $m2Chr){ $m1Type.="U";$m2Type.="D"; }
			
			}else{ 
				($UpRev,$DoRev)=($mateReverseStrand,$readReverseStrand);
				@ALen=($m2Len,$m1Len);
				@ASt=(($m2St-$cMin),($m1St-$cMin));
				
				if($m1Chr eq $m2Chr){ $m1Type.="D";$m2Type.="U"; }
			}
				
			if($m1Chr ne $m2Chr){ #orphan
				$mType="BND";
			}elsif($UpRev==0 and $DoRev==1){
				$mType="DEL"; # del, red
			}elsif($UpRev==1 and $DoRev==0){
				$mType="DUP"; # dup, blue
			}elsif($UpRev==$DoRev){ # fw fw and rw rw
				$mType="INV"; # inv, blue
			}else{die "s1 $_" }
		
			
			
			
			#### what are the breakpoint coordinates, make the breakpoint side appear larger in IGV   ---#   #---- or #---   ---#, # thickEnd only accepts one coord, so only possible for DEL and non connected features
			if($m1Ori eq "+"){   $m1Brkp=$m1En; $m1NBrkp=$m1St;  $m1ThickStEn=($m1En-5-1)."\t".$m1En;  }
			else{                $m1Brkp=$m1St; $m1NBrkp=$m1En;  $m1ThickStEn=($m1St-1)."\t".($m1St+5);  }
			if($m2Ori eq "+"){   $m2Brkp=$m2En; $m2NBrkp=$m2St;  $m2ThickStEn=($m2En-5-1)."\t".$m2En;  }
			else{                $m2Brkp=$m2St; $m2NBrkp=$m2En;  $m2ThickStEn=($m2St-1)."\t".($m2St+5);  }
	
				
			# sort breakpoint location now and at retrieval the same way to prevent the need to search for each breakpoint in case on different chromosomes
			@rp = sort { $$a[0] cmp $$b[0] || $$a[1] <=> $$b[1];  } ([$m1Chr,$m1Brkp,$m1Type,$m1NBrkp,$m1MQ],[$m2Chr,$m2Brkp,$m2Type,$m2NBrkp,$m2MQ]);
			
			$BRKPout= join("\t",($rp[0][0],$rp[0][1],"type=".join(",",( sort ($rp[0][2],$rp[1][2]) )).";type1=$rp[0][2];readID=$cID;NBP1=$rp[0][3];MQ1=$rp[0][4];chr2=$rp[1][0];BP2=$rp[1][1];type2=$rp[1][2];NBP2=$rp[1][3];MQ2=$rp[1][4]" ) )."\n";				
			
			$IGVout="";
			if( abs($s[8])<$maxLinkSize and $m1Chr eq $m2Chr){ # print connected
						
				$gffTags="Name2=".$mType."_$cGlobLen\_$m1Type,$m2Type;readID=$cID;distance=$cGlobLen;MQ1=$m1MQ;MQ2=$m2MQ;loc1=$m1Chr:$m1St:$m1Ori:$m1Cigar;loc2=$m2Chr:$m2St:$m2Ori:$m2Cigar;";	
				$IGVout.= join("\t",($m1Chr,($AMM[0]-1),$AMM[$#AMM],$gffTags,($m1MQ+$m2MQ/2),$m1Ori,($AMM[0]-1),$AMM[$#AMM],type2colFun("PE",$mType,$m2Chr),2,join(",",@ALen),join(",",@ASt)))."\n";
			
			}else{ # print pairs seperate
						
				$gffTags="Name2=".$mType."_$m1Type;readID=$cID;MQ=$m1MQ;loc2=$m2Chr:$m2St;cigar=$m1Ori:$m1Cigar;joinType=$mType;";
				$IGVout.= join("\t",($m1Chr,($m1St-1),$m1En,$gffTags,$m1MQ,$m1Ori,$m1ThickStEn,type2colFun("PE",$mType,$m2Chr),1,$m1Len,0))."\n";
			
				$gffTags="Name2=".$mType."_$m2Type;readID=$cID;MQ=$m2MQ;loc2=$m1Chr:$m1St;cigar=$m2Ori:$m2Cigar;joinType=$mType;";
				$IGVout.= join("\t",($m2Chr,($m2St-1),$m2En,$gffTags,$m2MQ,$m2Ori,$m2ThickStEn,type2colFun("PE",$mType,$m1Chr),1,$m2Len,0))."\n";				
			}
			
			
			push @aCluster, [$cID,$m1Chr,$m1Brkp,$m2Chr,$m2Brkp,0,1,$IGVout,$BRKPout,$_];
			clusterAndCleanup(1); ##### cluster on the fly, only output cluster of size >=2,  keep only 100 cluster at once
		
		}
	# 	last if $m1Chr ne "1";
	# 	last if ($c++)>100000;
	}
	
	clusterAndCleanup(0);
	
	close(IN1);
	
	if(substr($cMode,0,1) eq "d" ){
		print STDERR "number valid IDs ".scalar(keys %validDPIDs)."\n";
		########## outoutp SAM using the IDs, so that both reads of a pair end up in the bam file
		open(IN1, "samtools view $inFile | ") || die "  $inFile nicht gefunden";
		while(<IN1>){	
			chomp; @s=split("\t",$_);
			if (exists($validDPIDs{readID($s[0])})){
				print SAM $_."\n";
				$stat{passing_alignments}++;
			}
		}close(IN1);
	}
}



close(SAM);


foreach (sort { $stat{$b} <=> $stat{$a} } keys %stat){
	print STDERR "$_: ".$stat{$_}."\n";
}


sub clusterAndCleanup{
	
	($mode)=@_;
	
	
	if($mode==1){
		
		($AcID,$Am1Chr,$Am1Brkp,$Am2Chr,$Am2Brkp,$AclustID)=@{$aCluster[$#aCluster]};
		
# 		print STDERR "$AcID,$Am1Chr,$Am1Brkp,$Am2Chr,$Am2Brkp,$AclustID\n";
		#loop through existing cluster
		for($i=0; $i<$#aCluster; $i++){
		
			($BcID,$Bm1Chr,$Bm1Brkp,$Bm2Chr,$Bm2Brkp,$BclustID)=@{$aCluster[$i]};
		
			next if ($Am1Chr ne $Bm1Chr or $Am2Chr ne $Bm2Chr); # chr have to match
			
			next if ($Am1Brkp > ($Bm1Brkp+$fl) or $Am1Brkp < ($Bm1Brkp-$fl)  ); # pos 1 has to match
			
			next if ($Am2Brkp > ($Bm2Brkp+$fl) or $Am2Brkp < ($Bm2Brkp-$fl)  ); # pos 2 has to match
			
			
			$AclustID=$BclustID;
			$aCluster[$i][6]++;
			$aCluster[$#aCluster][6]++;
			
# 			print STDERR "## $i ## clsuter $BcID,$Bm1Chr,$Bm1Brkp,$Bm2Chr,$Bm2Brkp,$BclustID,  $Am1Chr=$Bm1Chr,  $Am2Chr=$Bm2Chr,  br1Diff:".(abs($Bm1Brkp-$Am1Brkp)).", br2Diff:".(abs($Bm2Brkp-$Am2Brkp))." clCount:$aCluster[$i][6]:$aCluster[$#aCluster][6]  \n";
		
			
		}
		
		if($AclustID == 0){
			$lClustID++; $aCluster[$#aCluster][5]=$lClustID;
		}else{
			$aCluster[$#aCluster][5]=$AclustID;
		}
	
	}
			
	# remove most distant cluster	
	while(  scalar(@aCluster)>100 or ($mode==0 and scalar(@aCluster)>0)){		
		
		print IGV $aCluster[0][7] if $aCluster[0][6] >=2;# and !exists($allowedChr{$aCluster[0][1]});
		print BRKP $aCluster[0][8] if $aCluster[0][6] >=2;
		
		$validDPIDs{$aCluster[0][0]}++ if $aCluster[0][6] >=2;
		
# 		print STDERR "## print ## count: $aCluster[0][6] ##  $aCluster[0][8]\n" if $aCluster[0][6] >=2;
# 		print STDERR "## not print ## count: $aCluster[0][6] ##  $aCluster[0][8]\n" if $aCluster[0][6] <2;
		
		
		shift @aCluster;		
	}
	
			
}
	


close(IGV);close(BRKP);
close(SAM);

print STDERR `sort_bgzip $outStem/bed/$cSample.$cMode.bed`;
print STDERR `sort_bgzip $outStem/bed/$cSample.$cMode.brkp`;


sub getSRType{

	my ($m1Chr,$m1MLoc,$m1Ori,$m1Dist,$m2Chr,$m2MLoc,$m2Ori,$m2Dist)=@_;
	
	my ($m1Type,$m2Type,$mType)=("","","");
	
	if($m1MLoc>0){
		($m1Type,$m2Type)=("1","2");
	}else{
		($m1Type,$m2Type)=("2","1");
	}

	$m1Type.=($m1Ori eq "+")? "p":"m";
	$m2Type.=($m2Ori eq "+")? "p":"m";
	
	if ($m1Chr eq $m2Chr){
		if($m1Dist>0){
			$m1Type.="U";
			$m2Type.="D";
			$mType="$m1Type,$m2Type";
		}else{
			$m1Type.="D";
			$m2Type.="U";
			$mType="$m2Type,$m1Type";	
		}
	}else{
		$mType=($m1MLoc>0)? "$m1Type,$m2Type":"$m2Type,$m1Type";
	}
	
	return ((($m1Type,$m2Type,$mType)));
		
}


sub getMatchLocationOnRead{
    my ($mOri,$mCigar)=@_;
	if($mOri eq "+"){
		if ($mCigar =~ /([0-9]+)[SH]$/){  return ($readLength-$1);   } # matching part is at the beginning of read 
		elsif($mCigar =~ /^([0-9]+)[SH]/){   return (($readLength-$1)*(-1));   } # matching part is at the end of read 
		else{ return 0}#die "s2 $mOri $mCigar $_"}
	}elsif($mOri eq "-"){
		if ($mCigar =~ /([0-9]+)[SH]$/){  return (($readLength-$1)*(-1))   } # matching part is at the beginning of read 
		elsif($mCigar =~ /^([0-9]+)[SH]/){  return ($readLength-$1);   } # matching part is at the end of read 
		else{ return 0}# die "s1 $mCigar"}
	}else{ die "s1 cigar $mCigar"}
}




print STDERR `tabix -f -p bed $outStem/bed/$cSample.$cMode.bed.gz`;
print STDERR `tabix -f -s 1 -b 2 -e 2  $outStem/bed/$cSample.$cMode.brkp.gz`;

sub readID{my ($rid)=@_; my @id=split(":",$rid); return join(":",@id[2..6])}
sub findEnd{my ($cigar)=@_; my $o=0; while ($cigar =~ m/([0-9]+)([MIDSH])/g){ ($L,$T)=($1,$2); if($T eq "M" or $T eq "D"){$o+=$L}} return $o}

sub clipPartCount{my ($cigar)=@_; my $o=0; while ($cigar =~ m/([0-9]+)([SH])/g){ $o++} return $o}



# 22	1000	5000	cloneA	960	+	1000	5000	195,67,67	2	567,488,	0,3512
# 22	2000	6000	cloneB	900	+	2000	6000	3,117,182	2	433,399,	0,3601
# 
# ST-E00185:49:H5LVWCCXX:5:1216:22264:23284_1     323     1       9999    3       107H43M 5       18606683        0       GATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC     KKKKKKKKKKKKKKKKKKKKKKKKKKKKFKKKKFFFKAKAAAK     NM:i:0  MD:Z:43 AS:i:43 XS:i:42 SA:Z:5,18606943,-,34S116M,60,0; XA:Z:hs37d5,-10061139,42M108S,0;        MC:Z:56M2D93M   MQ:i:60 RG:Z:H5LVWCCXX_5
# ST-E00185:49:H5LVWCCXX:5:1111:3253:71929_1      1395    1       9999    0       10H66M73H       5       18606943        0       GATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC      AKKAFKKKFFKKFFKAAKKKKAFKKKKKFKKKK7FKKKKKFKKK


