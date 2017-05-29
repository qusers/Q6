#!/bin/bash
#master script to check all energies of the basic tests
#runs one check for each possible combination of integrator and thermostat
#test is done for one of the three serial parallel hybrid
#tests if benchmark energies are present and screams if it is not the case

#Author: paul.bauer@icm.uu.se
#Released under Beer-Ware License

if [ $# -ne 1 ] ; then
echo "Number of arguments has to be one for this script to work"
echo "Please supply what test should be run"
echo "Supported arguments are serial parallel or hybrid"
exit 1
fi
if [ "$1" != "serial" ] && [ "$1" != "parallel" ] && [ "$1" != "hybrid" ] 
then
echo "Did not understand the argument $1"
echo "Please supply one of the following: serial parallel or hybrid"
exit 2
fi

wd=`pwd`

# Useful stuff:
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

###############################################################################
# Checking that the runs have finished without error.                         #
###############################################################################

function check_existing() {

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

}

function check_termination() {

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

}

function check_initial_energy() {

echo -n "Checking initial energy consistency                         "

INIT_QSURR=`sed -n '/Q-atom energies at step      0/,/==/ p' eq1.log | grep "Q-surr"` 
if [ "$1" == "SPH" ] ; then
INIT_QSURR_BM="Q-surr. 1 1.0000      3.12    139.43"
else
INIT_QSURR_BM="Q-surr. 1 1.0000    -31.30    228.64"
fi
if [ "$INIT_QSURR" == "$INIT_QSURR_BM" ]
then echo -e "$OK"
else 
 echo -e "$FAILED"
 echo "Initial Q-surr energy is: $INIT_QSURR"
 echo "Should be:                $INIT_QSURR_BM"
fi

}

###############################################################################
# Testing restart consistency. Checking that the final total potential energy #
# in the preceeding run equals the pot. energy in step 0 of the subsequent    #
# run.                                                                        #
# (Need not always be true, but in this case it is.)                          #
###############################################################################

function check_restart_consistency() {

sed -n '/Energy summary at step      0/,/==/ p' dc{2..5}.log \
 | grep "SUM"    > step0.tmp
sed -n '/FINAL  Energy summary/,/==/ p' dc{1..4}.log \
 | grep "SUM"   > final.tmp

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

}

###############################################################################
# Checking all energies. Produces a plot for gnuplot. There is                #
# an upper and lower bound from the benchmark runs of what the energies are   #
# expected to be. The further the simulation runs, the more random these      #
# energies become. The compilation should reproduce the energies in the       #
# first few steps and then fall somewhere in the ball park of the benchmark   #
# bounds. (Benchmark bounds were generated from a version of Q tested         #
# with the old benchmark by taking the average +/- std dev from twenty  runs  #
# with different random seeds. Hence this is a statistical figure of measure.)#
# you can make your own benchmark for a trusted version using the             #
# run_test.sh benchmark command It will run 20 different seeds by             #
# default and generate a new qsurr_benchmark.en file from them                #
###############################################################################

function check_all_energies() {

echo -n "Checking that benchmark energies are present                "

if [ -e $wd/${1}_${2}_${3}_benchmark.en ]
then echo -e "$OK"
else
 echo -e "$FAILED"
 echo "Could not locate benchmark energies for ${1} ${2} ${3}. Skipping..."
 exit 0
fi

echo -n "Extracting ligand-surrounding energies                      "

echo   '#          This Qdyn             Benchmark   '  > ${4}_${1}_${2}_${3}_NB.en
echo   '#          This Qdyn             Benchmark   '  > ${4}_${1}_${2}_${3}_QBND.en
echo   '#          This Qdyn             Benchmark   '  > ${4}_${1}_${2}_${3}_BND.en
echo   '#          This Qdyn             Benchmark   '  > ${4}_${1}_${2}_${3}_TOT.en
printf "# %10s %10s %10s %10s\n" "Qel" "Qvdw" "el" "vdw"  >> ${4}_${1}_${2}_${3}_NB.en
printf "# %10s %10s %10s %10s\n" "QBnd" "QAng" "QTor" "QImp" >> ${4}_${1}_${2}_${3}_QBND.en
printf "# %10s %10s %10s %10s\n" "Bnd" "Ang" "Tor" "Imp" >> ${4}_${1}_${2}_${3}_BND.en
printf "# %10s %10s %10s\n" "Tot" "Kin" "Pot" >> ${4}_${1}_${2}_${3}_TOT.en

grep "^Q-surr.*1 1.000" eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f \n",$4,$5}' >tmp1.en
grep "^solute " eq{1..5}.log dc{1..5}.log  | awk '{printf "%10.2f %10.2f \n",$2,$3}' >tmp2.en
grep "^solvent " eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f \n",$2,$3}' >tmp3.en
grep "^solute-solvent" eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f \n",$2,$3}' >tmp4.en
paste tmp2.en tmp3.en tmp4.en >tmp5.en
cat tmp5.en | awk '{printf "%10.2f %10.2f \n",$1+$3+$5,$2+$4+$6}' >tmp6.en
paste tmp1.en tmp6.en | nl >> ${4}_${1}_${2}_${3}_NB.en

/bin/rm tmp{1..6}.en

grep "^Q-atom " eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f %10.2f %10.2f \n",$4,$5,$6,$7}' | nl >> ${4}_${1}_${2}_${3}_QBND.en
grep "^solute " eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f %10.2f %10.2f \n",$4,$5,$6,$7}' | nl>> ${4}_${1}_${2}_${3}_BND.en
grep "^SUM " eq{1..5}.log dc{1..5}.log | awk '{printf "%10.2f %10.2f %10.2f \n",$2,$4,$3}' | nl >>${4}_${1}_${2}_${3}_TOT.en


echo -e "$OK"

}

function make_plots() {

echo -n "Preparing Gnuplot files 'NB_energy_benchmark.plot' 'QBND_energy_benchmark.plot' 'BND_energy_benchmark.plot' and 'TOT_energy_benchmark.plot'             "

scale="[kcal/mole]"

for i in NB QBND BND TOT ; do
if [[ "$i" == "TOT" ]] ; then
numplots=3
plotx=2
ploty=2
ylab1="Total energy $scale"
ylab2="Kinetic energy $scale"
ylab3="Potential energy $scale"
startn=24
else
numplots=4
plotx=2
ploty=2
if [[ "$i" == "NB" ]] ; then
ylab1="Q Nonbonded Electrostatics $scale"
ylab2="Q Nonbonded vdW $scale"
ylab3="Nonbonded Electrostatics $scale"
ylab4="Nonbonded vdW $scale"
startn=0
fi
if [[ "$i" == "QBND" ]] ; then
ylab1="Q Bonds $scale"
ylab2="Q Angles $scale"
ylab3="Q Torsions $scale"
ylab4="Q Impropers $scale"
startn=8
fi
if [[ "$i" == "BND" ]] ; then
ylab1="Bonds $scale"
ylab2="Angles $scale"
ylab3="Torsions $scale"
ylab4="Impropers $scale"
startn=16
fi
fi


#cat <<EOF > ${4}_${1}_${2}_${3}_${i}_qsurr.plot
echo "set title 'Plot for ${4} ${1} ${2} ${3} ${i}'
set multiplot layout ${plotx},${ploty}
set style line 1 lt 1 lc 1 lw 1 pt 6
set style line 2 lt 1 lc 3 lw 1 pt 6
set style line 3 lt 1 lc 5 lw 1 pt 6
"> ${4}_${1}_${2}_${3}_${i}_qsurr.plot
count=0
for ((jj=1;jj<=$numplots;jj++)) ; do
ylabel=ylab$jj
ylabel=$( echo $( eval echo \$$ylabel ))
onenum=$(( $jj + 1 + $startn + $count )) 
twonum=$(( $jj + 2 + $startn + $count ))

echo "
set xlabel 'Simulation step'
set ylabel '$ylabel'
set log x
set xrange [1:10000]
plot '$wd/${1}_${2}_${3}_benchmark.en' using 1:$onenum      \
                            notitle  with lines ls 1,       \
     '$wd/${1}_${2}_${3}_benchmark.en' using 1:$twonum      \
            title 'Benchmark bounds' with lines ls 1,       \
     '${4}_${1}_${2}_${3}_${i}.en' using 1:$(( ${jj} + 1 )) \
            title 'The present build' with lines ls 2
">>${4}_${1}_${2}_${3}_${i}_qsurr.plot
count=$(( count + 1 ))
done

echo -e "$OK"

echo -n "Checking if Gnuplot is installed on the system              "

if which gnuplot >& /dev/null 
then 
 echo -e "$OK"
 echo "Launching Gnuplot..."
 gnuplot -persist < ${4}_${1}_${2}_${3}_${i}_qsurr.plot &
else
 echo -e "$FAILED"
 echo "Gnuplot file qsurr.plot written."
fi

done 

}

#This part runs the actual analysis ...


for iii in SPH #  PBC 
do
for jjj in leap-frog # velocity-verlet
do
for ggg in berendsen # langevin nose-hoover
do
thisdir=`pwd`
if ! [ -f ${iii}_${jjj}_${ggg}_benchmark.en ] ; then
echo "No Benchmark energies present for ${iii} ${jjj} ${ggg}!"
echo "Run those first or get them from Github before running the tests"
continue
fi
cd run_${1}/${iii}/${jjj}/${ggg}
check_existing
check_termination
check_initial_energy $iii
check_restart_consistency
check_all_energies ${iii} ${jjj} ${ggg} $1
make_plots ${iii} ${jjj} ${ggg} $1
cd $thisdir
done
done
done

exit 0
