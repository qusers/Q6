#!/usr/bin/perl 
#fep_inper.pl
use warnings;
#       _____________________________________
#      |                                     |
#      |  created by Masoud Kazemi July 2011 |
#      |_____________________________________|
#	Modified by Paul Bauer March 2013
$a=0;
$pn=1;
while ($a<1.0001){ 
             $a= sprintf ( "%4.2f",$a);
             push (@landa, $a);
             $a=$a+0.02;
             $pn++;
}
$pn=100;
for ($b=0;$b<51;$b++){
$nn= sprintf ("%03d",$pn);
open INP, ">fep_$nn.inp";
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------#
print  INP "[MD]\n";
print  INP "steps                          50000\n";
print  INP "stepsize                        1.0\n";
print  INP "temperature                    300\n" ;
print  INP "separate_scaling                on\n" ;
print  INP "bath_coupling                     100\n";
print  INP "random_seed                     1250\n";
print  INP "lrf                               on\n";
print  INP "\n[cut-offs]\n";
print  INP "solute_solute                 10\n";
print  INP "solvent_solvent               10\n";
print  INP "q_atom                        99\n";
print  INP "\n[sphere]\n";                                       
print  INP "\n[water]\n";
print  INP "\n[intervals]\n"			      ;
print  INP "non_bond                      30\n"	      ;  # time step for nonbonded decreased
print  INP "output                        50\n"	      ;
print  INP "energy                         5\n"	      ;
print  INP "trajectory                   50\n"	      ;
print  INP "\n[files]\n"			      ;
print  INP "topology             co_0p50.top\n"        ;
printf INP "restart	         fep_%03d.re\n",$nn+2 ;
print  INP "final	         fep_$nn.re\n"	      ;
print  INP "energy	         fep_$nn.en\n"	      ;
print  INP "trajectory           fep_$nn.dcd\n"	      ;
print  INP "fep	                     co_0p50.fep\n"	      ;
print  INP "\n[lambdas]\n"			      ;
printf INP "%4.2f   %4.2f\n", (1-$landa[$b]) ,$landa[$b]   ;


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------#
$pn=$pn-2;
close INP ;
}
