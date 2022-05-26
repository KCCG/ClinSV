package My::Seq;
use strict;
use warnings;
 
use Exporter qw(import);
our @EXPORT_OK = qw(GC_seq);


sub GC_seq{


 my $str_pointer = shift;
 my $G=0;
 my $C=0;
 my $A=0;
 my $T=0;
 my $N=0;
 
 $G = $$str_pointer =~ tr/G//; # count the number of commas. this does not modify the string in any way
 $C = $$str_pointer =~ tr/C//;
 $A = $$str_pointer =~ tr/A//;
 $T = $$str_pointer =~ tr/T//;
 $N = $$str_pointer =~ tr/N//;
 
 my $AT=($A+$T);
 my $GC=($G+$C);
 
 
 if(($GC+$AT)==0){
 	return (-1);
 }
 
 return sprintf("%.3f",($GC/($AT+$GC)*100));
}








1;

