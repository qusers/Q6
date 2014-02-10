#!/usr/bin/perl
use strict;
use warnings;

open OUT , ">qG1P.prm" or die $!;

select OUT;

print "*------------------------------------------------\n";
print "*\n";
print "*Q-FF parameters: CHARMM27 parameters\n";
print "*------------------------------------------------\n";
print "[options] force-field options\n";
print "name            Q-charmm27 nucleic acids\n";
print "vdw_rule    arithmetic   ! vdW combination rule (geometric or arithmetic )\n";
print "scale_14         1.0         ! electrostatic 1-4 scaling factor\n";
print "switch_atoms  off	! on = use switch atom; off = use charge group\n";
print "improper_definition explicit ! improper representation by 2 or four atoms\n";
print "improper_potential	harmonic\n";
print "coulomb_constant 332.0716    ! Constant in electrostatic energy calculation; default = 332.\n";
print "force_field CHARMM           ! Force Field Type (GROMOS (default), AMBER or CHARMM)\n\n";
print "[atom_types] atom type definitions\n";
print "*-iac------Avdw1------Avdw2-----Bvdw1------Avdw3-----Bvdw2&3----mass---SYBYL-name-old-comment\n";

open IN1 , "<G1P.rtf" or die $!;
#open(IN1,$ARGV[0])  || die "Cannot open file \"$ARGV[0]\"";

my $counter1 = 0;
my @atomtype;
my @mass;

while (<IN1>) {
     if ($_ =~ /(MASS)\s+[0-9]+\s+(\w+)\s+([0-9]+\.[0-9]+)/) {
#     if ($_ =~ /(MASS)\s+0\.0\s+(\-[0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)/) {	 
          $atomtype[$counter1] = $2;
          $mass[$counter1] = $3;
          $counter1++;
     }
}

open IN2 , "<G1P.par" or die $!;
#open IN2 , "<>" or die $!;
#open(IN2,$ARGV[1])  || die "Cannot open file \"$ARGV[1]\"";

while (<IN2>) {
      if ($_ =~ /^(\w+)\s+0\.0\s+(\-[0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)/) {
          for my $i (0..$#atomtype) {
                 if ("$atomtype[$i]" eq "$1") {
                     printf "%-10s %-10.4f %-9.4f %-10.4f %-9.4f %-9.4f %6.3f %1.0f\n" , "$1" , $3 , $3 , (-1*$2) , $3 , (-1*$2) , $mass[$i] , 1;
                 }
          }
      }
}

print "*--------------------------------------------------\n\n";
print "[bonds] types definitions\n";
print "*iaci iacj   force.c.  dist.\n";
print "*------------------------------------------------\n";

open IN2 , "<G1P.par" or die $!;
#open IN2 , "<>" or die $!;
#open(IN2,$ARGV[1])  || die "Cannot open file \"$ARGV[1]\"";

while (<IN2>) {
      if ($_ =~ /^(\w+)\s+(\w+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)/) {
          printf "%-7s %-9s %8.3f %7.3f\n" , "$1" , "$2" , (2*$3) , $4;
      }
}

print "*------------------------------------------------\n\n";
print "[angles] angle type definitions\n";
print "* iaci iacj iack  forceK angle0\n";
print "*------------------------------------------------\n";

open IN2 , "<G1P.par" or die $!;
#open IN2 , "<>" or die $!;
#open(IN2,$ARGV[1])  || die "Cannot open file \"$ARGV[1]\"";

while (<IN2>) {
      if ($_ =~ /^(\w+)\s+(\w+)\s+(\w+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s*([0-9]+\.[0-9]+)?\s*([0-9]+\.[0-9]+)?/) {
          if ($6) {
             printf "%-7s %-7s %-9s %6.3f %7s %6s %10s\n" , "$1" , "$2" , "$3" , (2*$4) , "$5" , "$6" , "$7";
          }
          else {
             printf "%-7s %-7s %-9s %6.3f %7s %6s %10s\n" , "$1" , "$2" , "$3" , (2*$4) , "$5" , "0.00" , "0.00000";
          }
      }
}

print "*------------------------------------------------\n\n";
print "[torsions] torsion type definitions\n";
print "*iaci iacj iack iacl  forceK  #minima phase   #paths\n";
print "*------------------------------------------------\n";

#open IN2 , "<>" or die $!;
open IN2 , "<G1P.par" or die $!;
#open(IN2,$ARGV[1])  || die "Cannot open file \"$ARGV[1]\"";

while (<IN2>) {
     if ($_ =~ /^(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([0-9]+\.[0-9]+)\s+([0-9]+)\s+(\-?[0-9]+\.[0-9]+)/) {
          last if "$2" eq "X";
          if ("$1" eq "X") {
              printf "%-7s %-7s %-7s %-10s %-5.3f %7.3f %8.3f %5.0f\n" , "?" , "$2" , "$3" , "?" , $5 , $6 , $7 , 1;
          }
          else {
              printf "%-7s %-7s %-7s %-10s %-5.3f %7.3f %8.3f %5.0f\n" , "$1" , "$2" , "$3" , "$4" , $5 , $6 , $7 , 1;
          }
     }
}

print "*------------------------------------------------\n\n";
print "[impropers] improper torsion type definitions\n";
print "*iaci iacj forceK   imp0\n";
print "*------------------------------------------------\n";

#open IN2 , "<>" or die $!;
#open(IN2,$ARGV[1])  || die "Cannot open file \"$ARGV[1]\"";
open IN2 , "<G1P.par" or die $!;

while (<IN2>) {
    if ($_ =~ /(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\-?[0-9]+\.[0-9]+)\s+(0)\s+([0-9]+\.[0-9]+)/) {
         if ("$2" eq "X") {
              printf "%-7s %-7s %-7s %-8s %7.3f %8.3f\n" , "?" , "$1" , "$4" , "?" , (2*$5) , $7;
         }
         else {
              printf "%-7s %-7s %-7s %-8s %7.3f %8.3f\n" , "$1" , "$2" , "$3" , "$4" , (2*$5) , $7;
         }
    }
}
print "*------------------------------------------------\n";
