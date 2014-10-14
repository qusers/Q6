#!/bin/bash

# Useful stuff:
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

###############################################################################
# Checking that the runs have finished without error.                         #
###############################################################################

for step in {1..5}
do
 echo -n "Checking presence of equilibration step $step                   "

 if [ -e eq${step}.log ]
 then echo -e "$OK"
 else 
   echo -e "$FAILED"
   echo "eq${step}.log not present, run terminated? Aborting test..."
   exit 1
 fi
done

for step in {1..5}
do
 echo -n "Checking presence of production step $step                      "

 if [ -e dc${step}.log ]
 then echo -e "$OK"
 else 
   echo -e "$FAILED"
   echo "dc${step}.log not present, run terminated? Aborting test..."
   exit 1
 fi
done

for step in {1..5}
do
 echo -n "Checking normal termination in equilibration step $step         "

 if grep 'terminated normally' eq${step}.log >& /dev/null
 then echo -e "$OK"
 else 
   echo -e "$FAILED"
   echo "eq${step}.log did not terminate normally. Aborting test..."
   exit 1
 fi
done

for step in {1..5}
do
 echo -n "Checking normal termination in production step $step            "

 if grep 'terminated normally' dc${step}.log >& /dev/null
 then echo -e "$OK"
 else 
   echo -e "$FAILED"
   echo "dc${step}.log did not terminate normally. Aborting test..."
   exit 1
 fi
done

echo -n "Checking initial energy consistency                         "

INIT_QSURR=`sed -n '/Q-atom energies at step      0/,/==/ { /Q-surr./ p }' \
            eq1.log`
INIT_QSURR_BM="Q-surr. 1 1.0000      3.12    139.43"

if [ "$INIT_QSURR" == "$INIT_QSURR_BM" ]
then echo -e "$OK"
else 
 echo -e "$FAILED"
 echo "Initial Q-surr energy is: $INIT_QSURR"
 echo "Should be:                $INIT_QSURR_BM"
fi

###############################################################################
# Testing restart consistency. Checking that the final total potential energy #
# in the preceeding run equals the pot. energy in step 0 of the subsequent    #
# run.                                                                        #
# (Need not always be true, but in this case it is.)                          #
###############################################################################

sed -n '/Energy summary at step      0/,/==/ { /SUM/ p }' dc{2..5}.log \
    > step0.tmp
sed -n '/FINAL  Energy summary/,/==/ { /SUM/ p }' dc{1..4}.log \
    > final.tmp

step=1
paste -d' ' final.tmp step0.tmp | tr -s [:blank:] | \
while read row
do
 # Row example:
 # ---
 # SUM -8317.97 -10355.22 2037.26 SUM -8313.92 -10355.22 2041.30
 # SUM -8284.87 -10338.52 2053.65 SUM -8269.74 -10338.52 2068.78
 # ---          ^^^^^^^^^                      ^^^^^^^^^
 set $row
 previous_pot=$3
 next_pot=$7

 echo -n "Checking restart energy consistency, production step $step      "

 if [ "$previous_pot" == "$next_pot" ]
 then echo -e "$OK"
 else echo -e "$FAIL"
 fi

 (( step++ ))
done

rm {final,step0}.tmp

###############################################################################
# Checking ligand-surrounding energies. Produces a plot for gnuplot. There is #
# an upper and lower bound from the benchmark runs of what the energies are   #
# expected to be. The further the simulation runs, the more random these      #
# energies become. The compilation should reproduce the energies in the       #
# first few steps and then fall somewhere in the ball park of the benchmark   #
# bounds. (Benchmark bounds were generated from a 'trusted' version of        #
# Q by taking the average +/- std dev from ten runs with different random     #
# seeds. Hence this is a statistical figure of measure.)                      #
###############################################################################

echo -n "Checking that benchmark energies are present                "

if [ -e qsurr_benchmark.en ]
then echo -e "$OK"
else
 echo -e "$FAILED"
 echo "Could not locate benchmark energies. Skipping..."
 exit 0
fi

echo -n "Extracting ligand-surrounding energies                      "

echo   '#          This Qdyn             Benchmark   '  > this_qsurr.en
printf "# %10s %10s %10s %10s\n" "el" "vdw" "el" "vdw" >> this_qsurr.en
awk '/^Q-surr./ {printf "%10.2f %10.2f \n",$4,$5}' \
     eq{1..5}.log dc{1..5}.log | cat -n                >> this_qsurr.en

echo -e "$OK"

echo -n "Preparing Gnuplot file 'energy_benchmark.plot'              "

cat <<EOF > qsurr.plot
set multiplot layout 2,1
set style line 1 lt 1 lc 1 lw 1 pt 6
set style line 2 lt 1 lc 3 lw 1 pt 6
set style line 3 lt 1 lc 5 lw 1 pt 6

set xlabel 'Simulation step'
set ylabel 'VdW energy [kcal/mol]'
set log x
set xrange [1:3000]
plot 'qsurr_benchmark.en' using 1:4                         \
                            notitle  with lines ls 1,       \
     'qsurr_benchmark.en' using 1:5                         \
            title 'Benchmark bounds' with lines ls 1,       \
     'this_qsurr.en' using 1:3                              \
            title 'The present build' with lines ls 2

set ylabel 'El energy [kcal/mol]'
plot 'qsurr_benchmark.en' using 1:2                         \
                                 notitle  with lines ls 1,  \
     'qsurr_benchmark.en' using 1:3                         \
                 title 'Benchmark bounds' with lines ls 1,  \
     'this_qsurr.en' using 1:2                              \
          title 'The present build' with lines ls 2     
EOF

echo -e "$OK"

echo -n "Checking if Gnuplot is installed on the system              "

if which gnuplot >& /dev/null 
then 
 echo -e "$OK"
 echo "Launching Gnuplot..."
 gnuplot -persist < qsurr.plot &
else
 echo -e "$FAILED"
 echo "Gnuplot file qsurr.plot written."
fi

exit 0
