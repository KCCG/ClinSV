
$path=shift @ARGV;
$wigPATH="$path/wigs";



print STDERR "## concat the wig files\n";
##################################################

use File::Basename;
open(OUTS, ">$wigPATH/sample2pos.txt") || die "files[0] not found"; # open first file as an template

@files=glob("$wigPATH/*.wig");
for ($i=0; $i<=$#files;$i++){ # open the files
	open($fileH[$i], "<$files[$i]") || die "$files[$i] not found";
	print OUTS "$i\t".basename($files[$i])."\n";
}
open(OUTS);
open(OUT, ">$wigPATH/summary.tab") || die "files[0] not found"; # open first file as an template 


$chrPf='';
open(IN, "<$files[0]") || die "files[0] not found"; # open first file as an template
while(<IN>){
	if(/^f/){
		if(/^fixedStep chrom=([^ ]+) /){ $cChr=$1 }		
		else{die $_}
		print " $startPos $_ $cChr\n";
		$chrPf='chr' if $cChr eq "chr1";
		$startPos=1;
		@cFOrder= values %F;
		foreach $cFH (@fileH){ 
			$c=<$cFH>;
			if( $c =~ /^fixedStep chrom=([^ ]+) /){ $cChr2=$1 }
			else{ die "line $. no fixedStep found\n" }
			die "line $. $cChr2 ne $cChr \n" if $cChr2 ne $cChr;
		 } # skip the current line for all files
		next;
	}
	undef @t; foreach $cFH (@fileH){ $c=<$cFH>; chomp($c); push @t, $c;} # create teh current lines
	print OUT join("\t",($cChr,$startPos,@t))."\n";
	$startPos+=1000;		
}close(IN);
foreach(@files){ close($F{basename($_)});  }
close(OUT);

if ($chrPf=='chr'){
	$Y_chr="chrY";
	$X_chr="chrX";
	$M_chr="chrM";
	print STDERR "set chr prefix \n";
}else{
	$Y_chr="Y";
	$X_chr="X";
	$M_chr="MT";
}


print STDERR "## compute the average coverage\n";
##################################################
open(STSTOUT, ">$wigPATH/summary.cov_stat.tab") || die " $wigPATH/summary.tab not found";

open(IN1, "<$wigPATH/summary.tab") || die " $wigPATH/summary.tab not found";
while(<IN1>){ chomp; @_=split("\t",$_); 
	next if $_[0] eq "hs37d5" or $_[0] eq "NC_007605";
	
	if($_[0] =~ /^(chr)*[0-9]+$/){
		if ($_[1]>=20000001 and $_[1]<=30000000){
			for($i=2;$i<=$#_;$i++){
				$hS{A}{$i}+=$_[$i],$hC{A}{$i}++ if $_[$i]>0; 	
			}
		}
	}elsif($_[0] =~ /^(chr)*Y$/){
		if (($_[1]>=6641419 and $_[1]<=10079253 or $_[1]>=13800704 and $_[1]<=23668908)){ 
			for($i=2;$i<=$#_;$i++){
				$hS{$Y_chr}{$i}+=$_[$i],$hC{$Y_chr}{$i}++ if $_[$i]>0;    
			}
		}
	}else{
		for($i=2;$i<=$#_;$i++){
			$hS{$_[0]}{$i}+=$_[$i],$hC{$_[0]}{$i}++ if $_[$i]>0;
		} 
	}
	
}close(IN1);



for($i=2;$i<=$#_;$i++){ # calculate the mean
	
	foreach $cChr (keys %hS){
		if ($hC{$cChr}{$i}==0){
			print STDERR "$cChr $i 0\n";  $meanCov{$cChr}{$i}=(-1) 
		}else{
			$meanCov{$cChr}{$i}=$hS{$cChr}{$i}/$hC{$cChr}{$i};  
		}
	} # calculate mean
	
	## foreach $cChr (keys %hS){  if($_[0] !~ /^(A|X|Y|MT|[0-9]+)$/){  print STDERR "$_[0] "; $meanCov{$cChr}{$i}=$meanCov{A}{$i}   }  } # set all other non X Y MT to A coverage
	$meanXOri=round($meanCov{$X_chr}{$i},1);
	$meanYOri=round($meanCov{$Y_chr}{$i},1);
	$meanCov{$X_chr}{$i}=sprintf("%.0f",($meanCov{$X_chr}{$i}/$meanCov{A}{$i})/0.5)*0.5*$meanCov{A}{$i};
	$meanCov{$Y_chr}{$i}=0 if $meanCov{$Y_chr}{$i} == (-1); # if -1 then because it is not covered	
	$meanCov{$Y_chr}{$i}=sprintf("%.0f",($meanCov{$Y_chr}{$i}/$meanCov{A}{$i})/0.5)*0.5*$meanCov{A}{$i};
	$gender{$i}=($meanCov{$Y_chr}{$i}==0)? "F":"M"; push @male, $i if $meanCov{$Y_chr}{$i}!=0;
	print STSTOUT "average coverage sample\t$i\tgender\t$gender{$i}\tautosome\t".round($meanCov{A}{$i},1)."\tX\t".round($meanCov{$X_chr}{$i},1)."\tY\t".round($meanCov{$Y_chr}{$i},1)."\tMT\t".round($meanCov{MT}{$i},1)."\tX ori\t$meanXOri\tY ori\t$meanYOri\n";

}

close(STSTOUT);

print STDERR "## normalize\n";
##################################################

open(IN2, "<$wigPATH/summary.tab") || die " $wigPATH/summary.tab not found";
open(OUT, ">$wigPATH/summary.norm.tab") || die " $wigPATH/summary.tab not found";

while(<IN2>){ chomp; @_=split("\t",$_);
	if($_[0] =~ /^(chr)*[0-9]+$/){ $cChr="A"; }else{$cChr=$_[0]; }
	for($i=2;$i<=$#_;$i++){  if ($meanCov{$cChr}{$i}==0){ $_[$i]=(-2) }else{ $_[$i]=round($_[$i]/$meanCov{$cChr}{$i},2) if $_[$i]>0; } }
	if($_[0] eq $Y_chr){ print OUT join("\t",@_[(0,1,@male)])."\n"; }else{ print OUT join("\t",@_)."\n"; }
}
close(IN2);
close(OUT);


print STDERR "## split up by chromosome and shuffle order for each chromosome\n";
##################################################

use List::Util qw(shuffle);

open(IN1, "<$wigPATH/summary.norm.tab") || die " $wigPATH/summary.norm.tab not found";
while(<IN1>){ 
	chomp; @t=split("\t",$_); 
	
	if(!exists($fh{$t[0]})){
		open($fh{$t[0]},">$path/$t[0].tab");
		
		#create new order
		@cFOrder= (0,1,shuffle 2..$#t);
		print STDER "chr $t[0], new sample order: ".join("\t",@cFOrder)."\n";
	}
	
	print {$fh{$t[0]}} join("\t",map { $t[$_] } @cFOrder)."\n";
}

foreach $cFH (values %fh){ close($cFH) }



sub round { my $number = shift || 0; my $dec = 10 ** (shift || 0); return int( $dec * $number + .5 * ($number <=> 0)) / $dec;}