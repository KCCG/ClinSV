
# Author: Andre Minoche, a.minoche@garvan.org.au

# v2 adjust the cutoff of what a discordant insert size is, to keep the coverage low


# perl filterConcordantPairs-v2.pl /g/data2/gd7/research/Tony_SV-MA_20150506/alignments/FR07885193.bam H00MYCCXX_5 0.03:0.25 1:1000001-5000000

($inputBam,$covCutOff,$region,$refFasta)=@ARGV;

my($loX,$hiX)=split(/[:]/,$covCutOff);

# if ($loX > 1){ ($loX,$hiX)=(400,400); }

my($chr,$st,$en)=split(/[:\-]/,$region);
my $len=($en-$st+1);

@allInsertSizes=();



################# find insert size cutoff, so that the correct pairs does not surpass $adjustCoverageToX x
# the assumption is that correct and discordant pairs have a iSize transition zone

# find the discordant insert size cutoff based on all lanes in the bam
# in case different libraries with different fragment size selection was used, 
# it suboptimal for the small discordant inserts, so large deletions are not affected, 
# since I don't expect it to happen very often I won't make an exception for it



print STDERR "samtools view -T $refFasta $inputBam $region | \n";
open(IN1, "samtools view -T $refFasta $inputBam $region |") || die " nicht gefunden";
while(<IN1>){

	if(index($_,"@")==0){next;} @_=split("\t",$_); 
	#ST-E00106:140:H00MYCCXX:5:2219:30536:21333      163     1       10000   23      6S35M109S       =       10004   150     CTAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACGATAGGCAGAGCGCTAAGCCTGGCCATAGTGCTAAGCTTAGTGTAGATCTCGGCGGTAGCCGTATCACTAAAAAAAACACTAACCCCACAATTAATCATAACCCAAGAC  AAFFFKKKKKKKKKKKKKKKKKKKKKKKKKKKAKFKKKKK,,,F,,,7,,7,A,(,,,,,7,,(,(,,77,,7FKK7,,7FF,,,,,7,A7,,7((((<,AF,,A,<,,,,7,,,,,A7A,,<,,,,,
    
	next if ($_[1] & 8); # next if mate unmapped
	next if ($_[1] & 4); # next if read unmapped
	next if ! ($_[1] & 1); # next if read is not paired
	next if ($_[1] & 1024); # skip PCR duplicates
	next if ($_[1] & 256); # skip if not primary alignment
	
	
	next if $_[6] ne "="; # skip if reads of pair on different scaffolds
	next if $_[8]<0; # count each read pair once: only count the forward strand read
 
 	# the forward strand read must be upstream
	if($_[3]<$_[7]){ next unless !($_[1] & 16) and ($_[1] & 32)  } # if current read is upstream, then it must be the FW read an the other the reverse
	if($_[3]>$_[7]){ next unless ($_[1] & 16) and !($_[1] & 32)  } # if current read is the downstream read, it must be the reverse matching, and the other the forward matching


# 	print STDOUT $_[8]."\n";
	next if $_[8]>=1500;
	push @allInsertSizes, $_[8];

}
close(IN1);


$cSum=0; $hiC=0; foreach $_ (sort {$b <=> $a } @allInsertSizes){ $c1++; $hiC=$_; last if $c1 >= $hiX/1000000*$len }

$cSum=0; $loC=0; foreach $_ (sort {$a <=> $b } @allInsertSizes){ $c2++; $loC=$_; last if $c2 >= $loX/1000000*$len }

# print STDERR "insert cutoff low $loC, count $c2 (to obtain $loX x), high $hiC, count $c1 (to obtain $hiX x) test region length $len, region: $region, all iSizes:".scalar(@allInsertSizes)."\n";
print STDERR "insert cutoff low $loC (to obtain $loX rpm), high $hiC (to obtain $hiX rpm) test region length $len, region: $region\n";

exit 1 if $hiC==0;

open(IN1, "samtools view -T $refFasta -hF 3340 $inputBam | ") || die " nicht gefunden";
while(<IN1>){
if(index($_,"@")==0){print;next;} @_=split("\t",$_); 

next if $_[2] eq "hs37d5";
next if $_[2] eq "NC_007605";

if((!($_[1] & 16) and ($_[1] & 32) and $_[6] eq "=" and $_[3]<$_[7] and $loC<abs($_[8]) and abs($_[8])<$hiC ) or
(($_[1] & 16) and !($_[1] & 32) and $_[6] eq "=" and $_[3]>$_[7] and $loC<abs($_[8]) and abs($_[8])<$hiC )){next}

next if $_[5] !~ /[0-9]+M/;

print;
}
close(IN1);











