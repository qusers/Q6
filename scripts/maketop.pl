#!/usr/bin/perl
use strict; 
use warnings;

#if( @ARGV < 2) # is less than two arguments
#{
#    print "Usage: charmm36top_2_qtop.pl charmmtop qtop ";
#    exit 0;
#}
#else {
open IN , "<G1P.rtf" or die $!;
#open(IN,$ARGV[0])  || die "Cannot open file \"$ARGV[0]\"";
open OUT , ">G1P.qlib" or die $!;
#open(OUT,$ARGV[1]) || die "Cannot open file \"$ARGV[1]\"";
#}



select OUT;
print "*------------------------------------------------------------------\n";
print "*\n";
print "* Q residue library file: CHARMM36 library\n";
print "*\n";

my $activate_first = 1;
my $activate_second = 0;
my $activate_third = 0;
my $label = 1;
my $counter1 = 0;
my $counter2=1;
my $memoriser = 0;
my $res = "";
my @charge_groups;

while (<IN>) {
   $activate_first = 0 if $_ =~ /^(PRES)/;
   $activate_first = 1 if $_ =~ /\w+\s+(TIP3|SOD|MG)/;
   if ($_ =~ /^(RESI)\s+(\w+)/ and $activate_first) {
      print "*------------------------------------------------------------------\n";
      $label = 1;
      $activate_second = 1;
      $activate_third = 1;
      print "{$2}\n";
      printf "%3s %-6s\n" , "" , "[info]";
      printf "%7s %-10s %-11s %-6s %-12s %-4s\n" , "" , "SYBYLtype" , "RESIDUE" , "!SYBYL" , "substructure" , "type";
      printf "%3s %-6s\n" , "" , "[atoms]";
      $res = "$2";
   }
   if ($_ =~ /^(ATOM|GROUP)\s+(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\-?[0-9]+\.[0-9]+)?/ and $activate_first) {
      printf "%10.0f %-6s %-11s %6.3f\n", $label , "$2" , "$3" , $4 if $2 and $3 and $4;
      $charge_groups[$counter1] = $2 if $2;
      $charge_groups[$counter1] = $1 if "$1" eq "GROUP";
      $counter1++;
      $label++ unless "$1" eq "GROUP";
   }
   if ($_ =~ /^(BOND)\s+(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*(\w+\'?\'?)?\s*/ and $activate_first) {
      printf "%11s\n" , "[bonds]" if $activate_second == 1;
      $activate_second = 0;
      printf "%6s %-4s %-17s %15s\n" , "" , "$2" , "$3" , "!connect i -- j" if $2 and $3;
      printf "%6s %-4s %-17s %15s\n" , "" , "$4" , "$5" , "!connect i -- j" if $4 and $5;
      printf "%6s %-4s %-17s %15s\n" , "" , "$6" , "$7" , "!connect i -- j" if $6 and $7;
      printf "%6s %-4s %-17s %15s\n" , "" , "$8" , "$9" , "!connect i -- j" if $8 and $9;
      printf "%6s %-4s %-17s %15s\n" , "" , "$10" , "$11" , "!connect i -- j" if $10 and $11;
      $counter2++;
   }
   if ($_ =~ /^(IMPH)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s*(\w+)?\s*(\w+)?\s*(\w+)?\s*(\w+)?\s*(\w+)?\s*(\w+)?\s*(\w+)?\s*(\w+)?\s/ and $activate_first) {
       printf "%17s %35s\n" , "[connections]" , "!how to bond to previous and next" if $activate_third;
       printf "%12s %3s %-3s\n" , "head" , "" , "P" if $activate_third;
       printf "%12s %3s %-3s\n" , "tail" , "" , "O3'" if $activate_third;  
       printf "%15s\n" , "[impropers]" if $activate_third;  
       printf "%10s %6s %6s %6s\n" , "$2" , "$3" , "$4" , "$5";
       printf "%10s %6s %6s %6s\n" , "$6" , "$7" , "$8" , "$9" if $6 and $7 and $8 and $9;
       printf "%10s %6s %6s %6s\n" , "$10" , "$11" , "$12" , "$13" if $10 and $11 and $12 and $13;
       printf "%58s\n" , "[charge_groups] !charge groups, with switch atom first" if $activate_third and $memoriser !=75;
       if ($activate_third and $memoriser != 75) {
          for my $i ($memoriser..$#charge_groups) {
               print "    $charge_groups[$i]" unless "$charge_groups[$i]" eq "GROUP";
               print "\n" if "$charge_groups[$i]" eq "GROUP";                 
               last if "$charge_groups[$i]" eq "O3'";
          }
          $memoriser += $label+3 if $activate_third;
          print "\n";
       }
       if ($activate_third == 0 and $memoriser == 75) {
           printf "%58s\n" , "[charge_groups] !charge groups, with switch atom first";
           for my $i ($memoriser..$#charge_groups) {
               print "    $charge_groups[$i]" unless "$charge_groups[$i]" eq "GROUP";
               print "\n" if "$charge_groups[$i]" eq "GROUP";                 
               last if "$charge_groups[$i]" eq "O3'";
           }
          $memoriser += $label+3;
          print "\n";
       }
       $activate_third = 0;
   }
   if ("$_" =~ /PATCHING/) {
        printf "%58s\n" , "[charge_groups] !charge groups, with switch atom first" if "$res" eq "TIP3" or "$res" eq "SOD" or "$res" eq "MG";
        print "    OH2    H1    H2\n" if "$res" eq "TIP3";
        print "    SOD\n" if "$res" eq "SOD"; 
        print "    MG\n"  if "$res" eq "MG";
   }
}
print "*------------------------------------------------------------------\n";
