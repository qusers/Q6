#!/usr/bin/perl 
#fep_inper.pl
use warnings;
#       _____________________________________
#      |                                     |
#      |  created by Masoud Kazemi July 2011 |
#      |_____________________________________|
#       Modified by Paul Bauer 2013-2014
$a=1;
$pn=1;
while ($a>0.0001){ 
             $a= sprintf ( "%4.2f",$a);
             push (@landa, $a);
             $a=$a-0.02;
             $pn++;
}
$pn=100;
for ($b=0;$b<51;$b++){
$nn= sprintf ("%03d",$pn);
open INP, ">fep_$nn.inp";
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------#
print  INP "[MD]\n";
print  INP "steps                          5000\n";
print  INP "stepsize                        1.0\n";
print  INP "temperature                    300\n" ;
print  INP "bath_coupling                     100\n";
print  INP "separate_scaling                on\n" ;
print  INP "lrf                               on\n";
print  INP "\n[cut-offs]\n";
print  INP "solute_solute                 10\n";
print  INP "solvent_solvent               10\n";
print  INP "solute_solvent		   10\n";
print  INP "q_atom                        99\n";
print  INP "\n[sphere]\n";
print  INP "#excluded_freeze		on\n";                                       
print  INP "\n[water]\n";
print  INP "\n[intervals]\n"                              ;
print  INP "non_bond                      30\n"            ;  # time step for nonbonded decreased
print  INP "output                        100\n"          ;
print  INP "trajectory                    100\n"          ;
print  INP "energy                         10\n"            ;       ;
print  INP "!\n[group_contribution]\n";
print  INP "!residue all 153\n";
print  INP "!residue electro 153\n";
print  INP "!residue vdw 153\n";
print  INP "!residue full 234\n";
print  INP "!residue all 153 234\n";
print  INP "\n[files]\n"                          ;
print  INP "topology            2cjpFH_ionres_oplsa.top\n"        ;
printf INP "restart             fep_%03d.re\n",$nn+2 ;
print  INP "final               fep_$nn.re\n"             ;
print  INP "energy              fep_$nn.en\n"             ;
print  INP "trajectory          fep_$nn.dcd\n"             ;
print  INP "fep                 lig.fep\n"            ;
print  INP "restraint           fep_102_rest.re\n";
print  INP "\n[lambdas]\n"                        ;
printf INP "%4.2f   %4.2f\n", (1-$landa[$b]) ,$landa[$b]   ;
print  INP "\n[sequence_restraints]\n";
print  INP "5102   5124  0.5   0   0   0\n" ;
print  INP "\n[distance_restraints]\n";
print  INP "1629	5113	3.30	3.30	3.0	1\n" ;
print  INP "1629	5116	3.30	3.30	3.0	1\n" ;

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------#
$pn=$pn-2;
close INP ;
}
