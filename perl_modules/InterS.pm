

package My::InterS;

use List::BinarySearch::XS qw(binsearch_pos);

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(merge_segs);
 


############# initiate new object
# my $IntersObj=My::InterS->new;
sub new{
    my $class = shift;
    
    my $s = {
        Mtable => {},
        indexM => {},
        notSorted => 1,
        sCoord => ["",(-1),(-1),0,0.9,0], 
        iCoord => [(-1),(-1)]
    };
    # Print all the values just for clarification.
#     print STDERR "Mtable is $$s{Mtable}\n";
#     print STDERR "indexM is $$s{indexM}\n";
#     print STDERR "notSorted is $$s{notSorted}\n";
    bless $s, $class;
    return $s;
}

############# add to table

#$IntersObj->add($qChr,$qSt,$qEn);
sub add{
    my $s = shift;
    my $qChr = shift;
    my $qSt = shift;
    my $qEn = shift;
    my @annotL = @_;
    my $qLen = ($qEn-$qSt+1);
	($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	
# 	print STDERR "$qChr,$qSt,$qEn,$qLen, ".join(",",@annotL)."\n";
    push @{$$s{Mtable}{$qChr}}, [$qSt,$qEn,$qLen,@annotL];
    $$s{notSorted}=1;
}

############# query table
#@ovlLines=$IntersObj->queryL($qChr,$qSt,$qEn,0.9,100); #reciprocal,tolerance
sub query{
    my $s = shift;
    
    my ($qChr,$qSt,$qEn,$pcOvlCutoff,$ovlTolerance)=@_[0..4];
    ($qSt,$qEn)=($qEn,$qSt) if $qSt>$qEn;
	my $qLen = ($qEn-$qSt+1);
	
	@{$$s{sCoord}}=($qChr,$qSt,$qEn,$qLen,$pcOvlCutoff,$ovlTolerance);
    @{$$s{iCoord}}=((-1),(-1));
    
	my $Mtable=$$s{Mtable};
	my $indexM=$$s{indexM};
	
	if($$s{notSorted}){
# 		print STDERR "sort & index... $Mtable\n";
		
		# sort master table
		foreach my $aChr (values %$Mtable){ @$aChr=sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @$aChr; }
	
		# create index for master
		foreach my $xChr (keys %{$Mtable}){@{$$indexM{$xChr}}=(); my $iChr=$$indexM{$xChr}; foreach my $cL (@{$$Mtable{$xChr}}){ push @$iChr,$$cL[0]; } }
		
		$$s{notSorted}=0;
	}

	return 0 if !exists($$Mtable{$qChr});
	
# 	print STDERR "binsearch_range... ($qSt-$qLen),$qEn\n";
	@{$$s{iCoord}}=binsearch_range($$indexM{$qChr},($qSt-$qLen),$qEn);
# 	print STDERR "result: ".join("\t",@{$$s{iCoord}})."\n";	
}


############# next result line
sub next_line{

	my $s = shift;
	my ($qChr,$qSt,$qEn,$qLen,$pcOvlCutoff,$ovlTolerance)=@{$$s{sCoord}};
	my ($iSt,$iEn)=@{$$s{iCoord}};	
	my $retArr=();
	
	return 0 if $iSt==(-1) or $iSt>$iEn;
	my $aChr=$$s{Mtable}{$qChr}; # pointer to current chromosome array 
	
# 	print STDERR "$qSt,$qEn : $iSt,$iEn : ".scalar(@$aChr)."\n";
	
	for my $i ($iSt..$iEn){
		my $cL=$$aChr[$iSt]; # pointer to current line
		${$$s{iCoord}}[0]=$i+1;
		
		my ($tSt,$tEn,$tLen)=@$cL[0..2];
		my $cOvlLen=overlap($qSt,$qEn,$tSt,$tEn);
# 		print STDERR "$qSt,$qEn,$qLen  $tSt,$tEn,$tLen i:$i\n";

		if( (($cOvlLen/$qLen) >= $pcOvlCutoff or ($cOvlLen>0 and abs($cOvlLen-$qLen)<=$ovlTolerance))  
			and (($cOvlLen/$tLen) >= $pcOvlCutoff or ($cOvlLen>0 and abs($cOvlLen-$tLen)<=$ovlTolerance))  ){
	# 		if( ($cOvlLen/$qLen) >= $pcOvlCutoff  and ($cOvlLen/$tLen) >= $pcOvlCutoff ){
# 				print STDERR "ovl OK $iSt:  (($cOvlLen/$qLen) >= $pcOvlCutoff or ($cOvlLen>0 and abs($cOvlLen-$qLen)<=$ovlTolerance))  and (($cOvlLen/$tLen) >= $pcOvlCutoff or ($cOvlLen>0 and abs($cOvlLen-$tLen)<=$ovlTolerance))  \n";
				
				@$retArr=($cOvlLen,sprintf("%.2f",($cOvlLen/$qLen*100)),sprintf("%.2f",($cOvlLen/$tLen*100)),@$cL); #ovl, % cov on querry, % cov on target, tSt, tEn, tLen, tAnnot
				return $retArr;

		}
	}
	return 0;
	
	
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


sub merge_segs{ # Array must be: ($chr1, $st1, $en1, $x1, $y1, $z1 ... ), returned array is $chr, $st, $en, [$x1,$x2], [$y1,$y2], [$z1,$z2]
	my ($arr)=@_;
	my ($i, $j, $k)=(0,0,0);
		
	# convert to  $chr, $st, $en, [$x1], [$y1], [$z1], if not already done
	if (ref($$arr[0][3]) ne "ARRAY"){ 
		for ($i=0; $i<=$#$arr; $i++){
			for ($j=3; $j<=$#{$$arr[$i]};$j++){
				$$arr[$i][$j]=[$$arr[$i][$j]];
			}
		}
	}#else{print STDERR "array already converted\n"}
	
	my $intitialArraySize=scalar(@$arr);
	return $arr if $intitialArraySize <= 1;	
	
 	#do the merging
 	@$arr=sort {$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]} @$arr; # sort by chr, st, end

	my $merged=0;
	for ($i=0; $i<=$#$arr;$i++){	
		$j=$i+1;
		last if $j>$#$arr;
		next if $$arr[$i][0] ne $$arr[$j][0];
		
		if (overlap($$arr[$i][1],$$arr[$i][2],$$arr[$j][1],$$arr[$j][2])>0){
			
			my ($x1,$x2,$x3,$x4)=sort {$a <=> $b} ($$arr[$i][1],$$arr[$i][2],$$arr[$j][1],$$arr[$j][2]);
			$$arr[$i][1]=$x1;
			$$arr[$i][2]=$x4;
			for ($k=3; $k<=$#{$$arr[$j]};$k++){  push @{$$arr[$i][$k]}, @{$$arr[$j][$k]};  }
			splice(@$arr,$j,1);
			$i--; # since the array got spliced don't leave out the next field
			$merged++;
		}
# 		print "$#$arr\n";
	}
# 	print "merge array from $intitialArraySize to ".scalar(@$arr)."\n";
	merge_segs($arr) if $merged>0;
	return $merged;	
}


sub overlap{ 
	my($R1,$R2,$Q1,$Q2)=@_; 
	my ($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2); 
	my ($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);   
	my $ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; 
	my $returnVal=($ovlLen>0)? $ovlLen:0;return $returnVal;
}



1;





































