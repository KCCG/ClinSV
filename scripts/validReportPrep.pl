
# Author: Andre Minoche, a.minoche@garvan.org.au

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../perlib";  # use the parent directory
print STDERR "local perl lib: ".$FindBin::Bin."/../../perlib\n";

use My::VCF2;
use POSIX;
use Bio::DB::Big;
use File::Path qw(make_path remove_tree);
use File::Basename;

$projectDir=shift(@ARGV);
$inRef=shift(@ARGV);
$inRefStyle=shift(@ARGV);
$fileExt=shift(@ARGV);
$cSample=shift(@ARGV);
$nameStemJoinF=shift(@ARGV);

if ($inRefStyle =~ /chr/){
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

print STDERR "projectDir: $projectDir\n";
print STDERR "inRef: $inRef\n";
print STDERR "fileExt: $fileExt\n";
print STDERR "sample: $cSample\n";

#read in chr2len
open(IN1, "<$inRef.fai") || die "can not open $inRef";
while(<IN1>){ chomp; @_=split("\t",$_); $chr2len{$_[0]}=$_[1]; }close(IN1);


########################### define expected values, and range from NA12878

%stat;
@categories = ("bamfile","testRegion","PcNotPropPair","PcMapDiffChr","covStdev","concInsSizeMin","concInsSizeMax","insMean","insSD","inputDiscPairs","inputSplitReads","avgCovAutosomes","avgCovX","avgCovY","avgCovMT","avgYOrig");

@SvStatCats = ("svStatHigh","svStatPassHigh","svStatLow","svStatCNVs","svStatCNVLoss","svStatCNVGain","svStatNonCNV","svStatToolLumpy","svStatToolCNVnator","svStatToolBoth","svStatGene","svStatRare","svStatRarePassHigh","svStatRarePassHighGene","svStatRarePassHighGeneCNV");

@categories = (@categories,@SvStatCats); # add SV type categories

#################### per sample

	
	$doBamStats=1;
	
	if ($doBamStats){
		if(-e "$projectDir/alignments/$cSample/$cSample.bam" ){ $cBam="$projectDir/alignments/$cSample/$cSample.bam"; }
		elsif(-e "$projectDir/alignments/$cSample/$cSample.IR.bam" ){ $cBam="$projectDir/alignments/$cSample/$cSample.IR.bam"; }
		else{print STDERR "can not locate bam file\n"; exit 1; }
	
		print STDERR "#### current sample: $cSample \n";
		print STDERR "# use bam file: $cBam \n"; $stat{$cSample}{bamfile}=$cBam;
		############### from bam compute 
		# 
		# - coverage, STDEV
		# 
		# 
		# - % badmates

		# - insertsize modal, mean, stdev, plot
		$countDiffrChr=0;
		$countPcNotPropPair=0;
		$countAllPairs=0;
		undef @tmpAIns;
		make_path($projectDir."/SVs/qc/insertSizes") if (! -d $projectDir."/SVs/qc/insertSizes");
		open(OUTINS, ">$projectDir/SVs/qc/insertSizes/$cSample.txt") || die "can not read bam file";

		#$chrPf here refers to input bam chromsome naming convention
		open(INB, "samtools view -T $inRef $cBam ".$chrPf."1:20000001-30000000 | ") || die "can not read bam file"; 
		# ST-E00141:48:H0ATBALXX:1:1220:8105:62471        163     1       19999858        60      150M    =       20000040        332     TCCAGAGTCTGTCTTGT  <AFFFJJF  NM:i:1  MD:Z:149A0      AS:i:149        XS:i:70 MC:Z:150M       MQ:i:60 RG:Z:H0ATBALXX_1
		while(<INB>){
		
			chomp; @_=split("\t",$_);
		
			next if (!($_[1] & 64)); # next if not first in pair
		
			$countAllPairs++;
			
			if($_[1] & 2){  if($_[6] eq "="){ print OUTINS abs($_[8])."\n"; push @tmpAIns, abs($_[8])  if abs($_[8]) < 1500;  }   }    
		
			if(!($_[1] & 2)){  $countPcNotPropPair++; }
			if($_[6] ne "="){  $countDiffrChr++;  }
		


		
		}
		close(INB);
	
		close(OUTINS);
		print STDERR " of $countAllPairs read pairs in region ".$chrPf."1:20000001-30000000\n";
		$PcNotPropPair=round($countPcNotPropPair/$countAllPairs*100,3);
		$PcMapDiffChr=round($countDiffrChr/$countAllPairs*100,3);
		print STDERR " % not proper pair: $PcNotPropPair\n";
		print STDERR " % mapping different Chr: $PcMapDiffChr\n";

		$stat{$cSample}{testRegion}="".$chrPf."1:20000001-30000000";
		$stat{$cSample}{PcNotPropPair}=$PcNotPropPair;
		$stat{$cSample}{PcMapDiffChr}=$PcMapDiffChr;
	
		$stat{$cSample}{insMean}=round(average(\@tmpAIns),0);
		$stat{$cSample}{insSD}=round(stdev(tmpAIns),0);
	
	
	
		####### stdev
	
		undef @tmpArr;
		#This  refers to ref genome
		open(INB, "samtools depth --reference $inRef -r ".$chrPf."1:20000001-30000000 $cBam | ") || die "can not read bam file"; 
		# ST-E00141:48:H0ATBALXX:1:1220:8105:62471        163     1       19999858        60      150M    =       20000040        332     TCCAGAGTCTGTCTTGT  <AFFFJJF  NM:i:1  MD:Z:149A0      AS:i:149        XS:i:70 MC:Z:150M       MQ:i:60 RG:Z:H0ATBALXX_1
		while(<INB>){
			chomp; @_=split("\t",$_);
			push @tmpArr, $_[2];
		}
		close(INB);

		$covStdev=round(stdev(\@tmpArr),1);
		print STDERR " covStdev: $covStdev\n";
		$stat{$cSample}{covStdev}=$covStdev;
		undef @tmpArr;
	
	}
	
	############### from STDERR output 
	# 
	# - select used insS cutoffs

	if ( -e "$projectDir/SVs/$cSample/lumpy/sh/lumpy.preproc.$cSample.e"){
		$insSCutOffPERaw=`grep "insert cutoff low" $projectDir/SVs/$cSample/lumpy/sh/lumpy.preproc.$cSample.e`;
	}else{
		$insSCutOffPERaw=`grep "insert cutoff low" $projectDir/SVs/$cSample/lumpy/$cSample.lumpy.metrics`;
	}
	
	if($insSCutOffPERaw =~ /insert cutoff low ([0-9]+) \(to obtain [0-9]+ rpm\), high ([0-9]+)/ ){  $concPErange="$1-$2";   }
	print STDERR " used insertsize range for concordant pairs: $concPErange\n";
	($stat{$cSample}{concInsSizeMin},$stat{$cSample}{concInsSizeMax})=split(/\-/,$concPErange);

	############### SE,PR
	# 
	# - input split reads discordant pairs
	
	if(-e "$projectDir/SVs/$cSample/lumpy/$cSample.discordants.bam"){
		$inputPEs=`samtools view -c $projectDir/SVs/$cSample/lumpy/$cSample.discordants.bam`; chomp($inputPEs);
		print STDERR "# use discordants file $projectDir/SVs/$cSample/lumpy/$cSample.discordants.bam  \n";
	}else{ # in case of multiple lanes per sample
		foreach $cFile (glob("$projectDir/SVs/$cSample/lumpy/$cSample\_*.discordants.bam")){
			print STDERR "# use discordants file $cFile  \n";
			$inputPEsX=`samtools view -c $cFile`; chomp($inputPEsX);
			$inputPEs+=$inputPEsX;
		}
	}
	$inputSRs=`samtools view -c $projectDir/SVs/$cSample/lumpy/$cSample.splitters.f.bam`; chomp($inputSRs);
	
	print STDERR " number input discordant pairs: $inputPEs\n";
	print STDERR " number input split reads: $inputSRs\n";
	
	$stat{$cSample}{inputDiscPairs}=$inputPEs; 
	$stat{$cSample}{inputSplitReads}=$inputSRs; 


	############### chromosome wide plot
	# 
	# - 1Mb window, with STDEV from control samples as error bars 


	print STDERR "open the bw files...\n";
	############################################################
	our %bwObj;


	foreach $cT (("mq","q0","q20")){			
# 		$cQBW=$jobfs."/".$cSample.".$cT.bw";
		$cQBW=$projectDir."/alignments/$cSample/bw/".$cSample.".$cT.bw";
		die " $cQBW ! exists" if (! -f $cQBW);
		$bwObj{$cSample}{$cT} = Bio::DB::Big->open($cQBW);
	}	
	print STDERR "# sample: $cSample, bw found\n";


	print STDERR "determine the average coverage per sample acrosse >MQ50 regions...\n";
	############################################################
	%meanCov;

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

	$stat{$cSample}{avgCovAutosomes}=$meanCov{A}{$cSample}; 
	$stat{$cSample}{avgCovX}=$meanCov{$X_chr}{$cSample}; 
	$stat{$cSample}{avgCovY}=$meanCov{$Y_chr}{$cSample}; 
	$stat{$cSample}{avgCovMT}=$meanCov{$chrPf."MT"}{$cSample}; 
	$stat{$cSample}{avgYOrig}=$meanCovYOrig; 

		
	######### 1 dot per MB,  output chr:value,value,value
	make_path($projectDir."/SVs/qc/chromView") if (! -d $projectDir."/SVs/qc/chromView");
	open(OUTCHRCOV, ">$projectDir/SVs/qc/chromView/$cSample.txt") || die "can not read bam file";
	
	$num{$X_chr}=round(round($meanCov{$X_chr}{$cSample}/$meanCov{A}{$cSample},1)*2,0);
	$num{$Y_chr}=round(round($meanCov{$Y_chr}{$cSample}/$meanCov{A}{$cSample},1)*2,0);
	
	print STDERR "numX: ".$num{$X_chr}.", numY: ".$num{$Y_chr}.", ".("X"  x  $num{$X_chr} ).", ".("Y"  x  $num{$Y_chr} )."\n";

	
			
	for $cChrPre ((1..22),"X","Y","MT"){
		$cChr=$chrPf.$cChrPre;
		$cChr=$M_chr if $cChr eq "chrMT";
		
		$qBins=ceil($chr2len{$cChr}/1000000);
		undef @aCov; undef(@covVals);

		eval {
				my $test = $bwObj{$cSample}{"q0"}->get_stats($cChr, 0, $chr2len{$cChr}, $qBins, 'mean');
		} or do {
				print STDERR "$@\n";
				print STDERR "trying different chromosome naming convention: $cChrPre\n";
				$cChr = $cChrPre;
		};

		$cDivCov=$meanCov{$cChr}{$cSample};
		$cChr2=$cChr;
		
		if( $cChr eq $X_chr  or $cChr eq $Y_chr ){
			# new logic: if XX devide by autosomal average
			#            if X,Y0,X0    devide by 1/2 autosomal average
			$cDivCov=$meanCov{A}{$cSample}*0.5 if $num{$cChr}<2;	
			
			if ($num{$cChr}==1){
				$cChr2=$cChr
			}elsif ($num{$cChr}>1){
				print "xxxxxx $cChr ".(substr($cChr,-1) x ($num{$cChr}-1) );
				$cChr2=$cChr.(substr($cChr,-1) x ($num{$cChr}-1) );
			}else{ # ==0
				$cChr2=$cChr."0"; # XX or X or Y or X0 or Yo
			}
			
		}
		
		my $aCov_stat = $bwObj{$cSample}{"q0"}->get_stats($cChr, 0, $chr2len{$cChr}, $qBins, 'mean');

		for ($i=0; $i<=$#$aCov_stat;$i++){		
		
		  $cCov=${$aCov_stat}[$i];
		  
		  $cCovVals=round($cCov/$cDivCov,2);
		  $cCovVals+=1 if  ($cChr eq $X_chr  or $cChr eq $Y_chr ) and $num{$cChr} == 0 and $cCovVals>=0;
		  push @covVals, $cCovVals;
		}
		
		print OUTCHRCOV "$cChr2\t".join("\t",@covVals)."\n";
	}
	close(OUTCHRCOV);

	############### number variants
	# 
	# type: CNVs, LOSS, GAIN, non CNV
	# 
	# tool: CNVnator, Lumpy, merged calls
	# 
	# PASS:  LOW, PASS, HIGH
	# 
	# hitting genes
	# 
	# RARE variants
	# 
	# RARE PASS HIGH
	# 
	# RARE PASS HIGH hitting genes
	# 
	# list of these variants # with R
	# DEL, DUP, other, sorted by size or evidence
	# size distribution, expected in grey
	# 

	

	$inVCF="$projectDir/SVs/$nameStemJoinF/SV-CNV$fileExt.vcf";
	print STDERR "inVCF $inVCF\n";
	my $V=My::VCF2->new($inVCF);

	while (my $V=$V->next_line){
		
		next if exists($$V{INFO}{SECONDARY});
		
		$stat{$cSample}{svStatHigh}++ if $$V{GTS}{FT}{$cSample} eq "HIGH";
		$stat{$cSample}{svStatPassHigh}++ if $$V{GTS}{FT}{$cSample} eq "PASS" or $$V{GTS}{FT}{$cSample} eq "HIGH";
		$stat{$cSample}{svStatLow}++ if $$V{GTS}{FT}{$cSample} eq "LOW";
		
		next if $$V{GTS}{FT}{$cSample} ne "PASS" and $$V{GTS}{FT}{$cSample} ne "HIGH";
		
		if ( $$V{INFO}{CNV}==1 or ( exists($$V{INFO}{CNEUTR}) and $$V{INFO}{CNEUTR}==0 ) ){
			$stat{$cSample}{svStatCNVs}++;
			if ($$V{INFO}{SVTYPE} eq "DEL"){
				$stat{$cSample}{svStatCNVLoss}++;
			}elsif($$V{INFO}{SVTYPE} eq "DUP"){
				$stat{$cSample}{svStatCNVGain}++;
			}elsif($$V{GTS}{DRA}{$cSample} < 1){
				$stat{$cSample}{svStatCNVLoss}++;
			}else{
				$stat{$cSample}{svStatCNVGain}++;
			}
		}else{  $stat{$cSample}{svStatNonCNV}++; } # next; }
				
		$stat{$cSample}{svStatToolLumpy}++ if $$V{INFO}{TOOL} eq "Lumpy";
		$stat{$cSample}{svStatToolCNVnator}++ if $$V{INFO}{TOOL} eq "CNVnator";
		$stat{$cSample}{svStatToolBoth}++ if $$V{INFO}{TOOL} eq "Lumpy,CNVnator";			
			
		$stat{$cSample}{svStatGene}++ if $$V{INFO}{NUMG}>0;
				
		$stat{$cSample}{svStatRare}++ if $$V{GTS}{RARE}{$cSample} ==1;
		
		$stat{$cSample}{svStatRarePassHighGene}++ if $$V{GTS}{RARE}{$cSample} ==1 and  $$V{INFO}{NUMG}>0;
		$stat{$cSample}{svStatRarePassHighGeneCNV}++ if $$V{GTS}{RARE}{$cSample} ==1 and  $$V{INFO}{NUMG}>0 and ( $$V{INFO}{CNV}==1 or ( exists($$V{INFO}{CNEUTR}) and $$V{INFO}{CNEUTR}==0 ) );
		
		$stat{$cSample}{svStatRarePassHigh}++ if $$V{GTS}{RARE}{$cSample} ==1 ;

	}
	

	


	############### sensitivity with NA12878 gold standard
	# 




	############### reprod CNVs with NA12878
	# 
	
	# prepare the 







############### write the variables for each sample into a table to be read with R

open(OUTTAB, ">$projectDir/SVs/qc/$cSample.stat.tab") || die "can not read bam file";
print OUTTAB join("\t",("Categroy",($cSample)))."\n";
foreach $cCat (@categories){

	
	
	print OUTTAB join("\t",($cCat, map { (exists($stat{$_}{$cCat}))? $stat{$_}{$cCat}:0 } ($cSample) ))."\n";
}

if (basename($inRef) =~ /38/){
	print OUTTAB join("\t",('build','38'))."\n";
}else{
	print OUTTAB join("\t",('build','37'))."\n";
}
close(OUTTAB);


######################################################################  function


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

	my $aCov_stat = $bwObj{$cSample}{$covQ}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean');
	my $aMQ_stat = $bwObj{$cSample}{mq}->get_stats($qChr, $qSt, $qEn, $qBins, 'mean') if $qMQ>0;	
	
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



sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}




sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}




































