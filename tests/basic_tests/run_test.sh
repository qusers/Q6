#!/bin/bash -l
#Script to submit/run basic binary tests for qdyn
#Arguments are serial/parallel/hybrid/benchmark
#to run the respective tests for all integrators and
#thermostats that qdyn supports, both for PBC and SPH
#Each run is analysed using the eval_test.sh script in the created
#subdirectory to compare to the benchmark results
#The following will be plotted: Total system energy, KIN energy,
#POT energy, Q-NB energy, NB-energy
#Q-BND energy, Q-ANG energy, Q-TOR energy, Q-IMP energy
#BND energy, ANG energy, TOR energy, IMP energy


QPATH="/data/simulations/Q_metal/q_github_paul/bin"
MODULES="intel/14.0 openmpi/1.8.1"
CORES=4
OMPCORES=2
PROJECT="SNIC2014-11-2"
MAIL="paul.bauer@icm.uu.se"
BENCHMARKS=20

if [ $# -ne 1 ] ; then
echo "Number of arguments has to be one for this script to work"
echo "Please supply what test should be run"
echo "Supported arguments are serial parallel hybrid and benchmark"
exit 1
fi
if [ "$1" != "serial" ] && [ "$1" != "parallel" ] && [ "$1" != "hybrid" ] && [ "$1" != "benchmark" ]
then
echo "Did not understand the argument $1"
echo "Please supply one of the following: serial parallel hybrid and benchmark"
exit 2
fi

#test if sbatch is there or not for later submission of jobs
findsbatch=`which sbatch`
if [ -z $findsbatch ] ; then
subcommand="bash "
else
subcommand="sbatch "
fi

#start of function block

#writes the file header for submission of the job
#arguments are file name and qbin type
function submit_header() {
echo "#!/bin/bash -l
#SBATCH -J bintest
#SBATCH -n $CORES
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH -A $PROJECT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$MAIL

module load $MODULES
" >$1
}

#This makes the actual run test script part
#arguments are file name and type of run
function make_run_test() {
if [ "$2" == "serial" ] || [ "$2" == "benchmark" ] ; then
RUNCOME="\$qbinary eq\${step}.inp > eq\${step}.log"
RUNCOMD="\$qbinary dc\${step}.inp > dc\${step}.log"
QBIN=qdyn5
fi
if [ "$2" == "parallel" ] ; then
RUNCOME="$MPICOM -np \$CORES \$qbinary eq\${step}.inp > eq\${step}.log"
RUNCOMD="$MPICOM -np \$CORES \$qbinary dc\${step}.inp > dc\${step}.log"
QBIN=qdyn5p
fi
if [ "$2" == "hybrid" ] ; then
RUNCOME="$MPICOM -np 1 -x OMP_NUM_THREADS=2 \$qbinary eq\${step}.inp > eq\${step}.log : -np \$(( \$CORES - 2 )) -x OMP_NUM_THREADS=1 \$qbinary"
RUNCOMD="$MPICOM -np 1 -x OMP_NUM_THREADS=2 \$qbinary dc\${step}.inp > dc\${step}.log : -np \$(( \$CORES - 2 )) -x OMP_NUM_THREADS=1 \$qbinary"
QBIN=qdyn5h
fi

echo "
QDIR=$QPATH
qbinary=\$QDIR/$QBIN
wd=\`pwd\`

if [ -z \$qbinary ]
then
 echo \"Please set the qbinary variable to point to the Q folder\"
 exit 1
elif [ ! -x \$qbinary ]
then
 echo \"Can't locate $QBIN in the \$qbinary variable, or you don't have
       execute permisson.\"
 exit 1
else
 echo \"Detected $QBIN in \${QDIR}\"
fi
# Useful vars
OK=\"(\\033[0;32m   OK   \\033[0m)\"
FAILED=\"(\\033[0;31m FAILED \\033[0m)\"

if [[ \$SNIC_TMP != \"\" ]];then
MYTMPDIR=\$SNIC_TMP/run_\$\$_$QBIN


else
#Most likely a local run
MYTMPDIR=\$\$
MYTMPDIR=/dev/shm/\${MYTMPDIR}_run_$QBIN

fi

rm -f eq{1..5}.log dc{1..5}.log >& /dev/null
tdir=\"\${MYTMPDIR}\"
cdir=\`pwd\`
for step in {1..5}
do
 mkdir \${tdir}
 echo \"Running equilibration step \${step}\"
 cp eq\${step}.inp *fep *top *re \${tdir}
 cd \${tdir}
 echo \"Entering temporary directory\"
 echo \"Now in \${tdir}\"
 $RUNCOME
 check=\$( grep \"terminated normally\" eq\${step}.log | wc -l)
 if [ ! \$check -eq 0 ]
 then echo -e \"\$OK\\nCopying back the files\"
 cp * \${cdir}
 cd \${cdir}
 rm -r \${tdir}
 else
  echo -e \"\$FAILED\"
  echo \"Check output (eq\${step}.log) for more info.\"
  echo \"Copying back the files\"
  cp * \${cdir}
  cd \${cdir}
  rm -r \${tdir}
  exit 1
 fi
done

for step in {1..5}
do
 mkdir \${tdir}
 echo \"Running production run step \${step} of 5\"
 cp dc\${step}.inp *fep *top *re \${tdir}
 cd \${tdir}
 echo \"Entering temporary directory\"
 echo \"Now in \${tdir}\"
 $RUNCOMD
 check=\$( grep \"terminated normally\" dc\${step}.log | wc -l)
 if [ ! \$check -eq 0 ]
 then echo -e \"\$OK\\nCopying back the files\"
 cp * \${cdir}
 cd \${cdir}
 rm -r \${tdir}
 else
  echo -e \"\$FAILED\"
  echo \"Check output (dc\${step}.log) for more info\"
  echo \"Copying back the files\"
  cp * \${cdir}
  cd \${cdir}
  rm -r \${tdir}
  exit 1
 fi
done">>$1
}

#this function modifies the inputs for the different runs
#to run the different thermostats and integrators
#as well as PBC and SPH boundaries
#argument is the type of calculation, with the list given later
#second one is either pbc or sph
function modify_inputs() {
if [ "$3" == "PBC" ] ; then
sed -i s/"invalid-boundary"/"\[PBC\]\nput_solvent_back_in_box\ton\nput_solute_back_in_box\toff"/g *inp
sed -i s/"invalid-qatom"/"q_atom\t24"/g *inp
sed -i s/"invalid-lrf"/"lrf\t24"/g *inp
else
sed -i s/"invalid-boundary"/"\[sphere\]\nshell_radius\t16\nshell_force\t20"/g *inp
sed -i s/"invalid-qatom"/"q_atom\t99"/g *inp
sed -i s/"invalid-lrf"/"lrf\t99"/g *inp
fi

sed -i s/"invalid-integrator"/"$1"/g *inp
sed -i s/"invalid-thermostat"/"$2"/g *inp

}

#function that creates the different directories and places the files in them
#argument is the run type serial parallel hybrid benchmark
function make_dirs() {
for i in SPH PBC ; do
for j in leap-frog velocity-verlet ; do
for g in berendsen langevin nose-hoover ; do
thisdir=`pwd`
if [ "$1" != "benchmark" ] ; then
if ! [ -f ${i}_${j}_${g}_benchmark.en ] ; then
echo "No Benchmark energies present for ${i} ${j} ${g}!"
echo "Run those first or get them from Github before running the tests"
continue
fi
mkdir -p run_${1}/${i}/${j}/${g}
cp *inp prep_${i}/*top prep_${i}/*fep run_${1}/${i}/${j}/${g}
cd run_${1}/${i}/${j}/${g}
modify_inputs $j $g $i
submit_header basetest.sh $1
make_run_test basetest.sh $1
$subcommand basetest.sh
cd $thisdir
else
CORES=`grep processor /proc/cpuinfo | wc -l`
for ((hh=1;hh<=$BENCHMARKS;hh++)) ; do
mkdir -p run_${1}/${i}/${j}/${g}/${hh}_benchmark
cp *inp prep_${i}/*top prep_${i}/*fep run_${1}/${i}/${j}/${g}/${hh}_benchmark
cd run_${1}/${i}/${j}/${g}/${hh}_benchmark
modify_inputs $j $g $i
sed -i s/"random_seed.*1"/"random_seed                $RANDOM"/g eq1.inp
echo "#!/bin/bash -l">run_${1}_${hh}.sh
make_run_test run_${1}_${hh}.sh $1
sed -i s/"run_"/"run_${hh}_"/g run_${1}_${hh}.sh
cd $thisdir
done
cd run_${1}/${i}/${j}/${g}/
submit_header basetest.sh $1
echo "
RUNNING=0
PARALLEL=$CORES
for ((run=1;run<=$BENCHMARKS;run++)) ; do
cd \${run}_benchmark

bash run_${1}_\${run}.sh &>logfile &
RARRAY[\$RUNNING]=\$\$
RUNNING=\$(( \${RUNNING} + 1 ))
cd ..

if [ \$RUNNING -eq \$PARALLEL ] || [ \$run -eq $BENCHMARKS ] ; then
wait
for ((testr=\$run;testr>=\$((\$run - \$RUNNING + 1));testr--)) ; do
        err=\`grep \"terminated normally\" \${testr}_benchmark/dc5.log | wc -l\`
        if [ \$err -eq 0 ] ; then
        echo \"Error during running of \${testr}_benchmark\"
        echo \"Aborting tests\"
        exit 666
        fi
done

RUNNING=0
fi

done
#Now get energies in nice arrays to prepare the new benchmark file
LENGTH=\`grep \"Q-surr. 1 1.0000\" 1_benchmark/eq*log 1_benchmark/dc*log | sed s/\".*Q-surr. 1 1.0000\"//g | awk '{print \$1}' | wc -l\`
for ((ave=0;ave<\$LENGTH;ave++)) ; do
AVQEL[\$ave]=0
STQEL[\$ave]=0
AVQLJ[\$ave]=0
STQLJ[\$ave]=0
AVEL[\$ave]=0
STEL[\$ave]=0
AVLJ[\$ave]=0
STLJ[\$ave]=0
AVQBND[\$ave]=0
STQBND[\$ave]=0
AVQANG[\$ave]=0
STQANG[\$ave]=0
AVQTOR[\$ave]=0
STQTOR[\$ave]=0
AVQIMP[\$ave]=0
STQIMP[\$ave]=0
AVBND[\$ave]=0
STBND[\$ave]=0
AVANG[\$ave]=0
STANG[\$ave]=0
AVTOR[\$ave]=0
STTOR[\$ave]=0
AVIMP[\$ave]=0
STIMP[\$ave]=0
AVTOT[\$ave]=0
STTOT[\$ave]=0
AVKIN[\$ave]=0
STKIN[\$ave]=0
AVPOT[\$ave]=0
STPOT[\$ave]=0
done

for ((run=1;run<=$BENCHMARKS;run++)) ; do
TARRAYQEL=( \`grep \"Q-surr.*1 1.0000\" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\".*Q-surr. 1 1.0000\"//g | awk '{print \$1}'\` )
TARRAYQLJ=( \`grep \"Q-surr.*1 1.0000\" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\".*Q-surr. 1 1.0000\"//g | awk '{print \$2}'\` )
TARRAYEL1=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$1}'\`)
TARRAYEL2=( \`grep \"^solvent \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solvent\"//g | awk '{print \$1}'\`)
TARRAYEL3=( \`grep \"^solute-solvent \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute-solvent\"//g | awk '{print \$1}'\`)
TARRAYLJ1=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$2}'\`)
TARRAYLJ2=( \`grep \"^solvent \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solvent\"//g | awk '{print \$2}'\`)
TARRAYLJ3=( \`grep \"^solute-solvent \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute-solvent\"//g | awk '{print \$2}'\`)
TARRAYQBND=( \`grep \"^Q-atom \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"Q-atom\"//g | awk '{print \$3}'\`)
TARRAYQANG=( \`grep \"^Q-atom \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"Q-atom\"//g | awk '{print \$4}'\`)
TARRAYQTOR=( \`grep \"^Q-atom \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"Q-atom\"//g | awk '{print \$5}'\`)
TARRAYQIMP=( \`grep \"^Q-atom \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"Q-atom\"//g | awk '{print \$6}'\`)
TARRAYBND=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$3}'\`)
TARRAYANG=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$4}'\`)
TARRAYTOR=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$5}'\`)
TARRAYIMP=( \`grep \"^solute \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log | sed s/\"solute\"//g | awk '{print \$6}'\`)
TARRAYTOT=( \`grep \"^SUM \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log  | awk '{print \$2}'\`)
TARRAYKIN=( \`grep \"^SUM \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log  | awk '{print \$4}'\`)
TARRAYPOT=( \`grep \"^SUM \" \${run}_benchmark/eq*log \${run}_benchmark/dc*log  | awk '{print \$3}'\`)


for ((ave=0;ave<\$LENGTH;ave++)) ; do
EARRAYQEL[\$run,\$ave]=\${TARRAYQEL[\$ave]}
EARRAYQLJ[\$run,\$ave]=\${TARRAYQLJ[\$ave]}
EARRAYEL[\$run,\$ave]=\`echo \${TARRAYEL1[\$ave]} \${TARRAYEL2[\$ave]} \${TARRAYEL3[\$ave]} | awk '{print \$1+\$2+\$3}'\`
EARRAYLJ[\$run,\$ave]=\`echo \${TARRAYLJ1[\$ave]} \${TARRAYLJ2[\$ave]} \${TARRAYLJ3[\$ave]} | awk '{print \$1+\$2+\$3}'\`
EARRAYQBND[\$run,\$ave]=\${TARRAYQBND[\$ave]}
EARRAYQANG[\$run,\$ave]=\${TARRAYQANG[\$ave]}
EARRAYQTOR[\$run,\$ave]=\${TARRAYQTOR[\$ave]}
EARRAYQIMP[\$run,\$ave]=\${TARRAYQIMP[\$ave]}
EARRAYBND[\$run,\$ave]=\${TARRAYBND[\$ave]}
EARRAYANG[\$run,\$ave]=\${TARRAYANG[\$ave]}
EARRAYTOR[\$run,\$ave]=\${TARRAYTOR[\$ave]}
EARRAYIMP[\$run,\$ave]=\${TARRAYIMP[\$ave]}
EARRAYTOT[\$run,\$ave]=\${TARRAYTOR[\$ave]}
EARRAYKIN[\$run,\$ave]=\${TARRAYKIN[\$ave]}
EARRAYPOT[\$run,\$ave]=\${TARRAYPOT[\$ave]}
done
for ((ave=0;ave<\$LENGTH;ave++)) ; do
AVQEL[\$ave]=\$( echo \"scale=5;\${AVQEL[\$ave]} + \${EARRAYQEL[\$run,\$ave]} \" | bc -q)
STQEL[\$ave]=\$( echo \"scale=5;\${STQEL[\$ave]} + (\${EARRAYQEL[\$run,\$ave]}*\${EARRAYQEL[\$run,\$ave]}) \" | bc -q)
AVQLJ[\$ave]=\$( echo \"scale=5;\${AVQLJ[\$ave]} + \${EARRAYQLJ[\$run,\$ave]} \" | bc -q)
STQLJ[\$ave]=\$( echo \"scale=5;\${STQLJ[\$ave]} + (\${EARRAYQLJ[\$run,\$ave]}*\${EARRAYQLJ[\$run,\$ave]}) \" | bc -q)
AVEL[\$ave]=\$( echo \"scale=5;\${AVEL[\$ave]} + \${EARRAYEL[\$run,\$ave]} \" | bc -q)
STEL[\$ave]=\$( echo \"scale=5;\${STEL[\$ave]} + (\${EARRAYEL[\$run,\$ave]}*\${EARRAYEL[\$run,\$ave]}) \" | bc -q)
AVLJ[\$ave]=\$( echo \"scale=5;\${AVLJ[\$ave]} + \${EARRAYLJ[\$run,\$ave]} \" | bc -q)
STLJ[\$ave]=\$( echo \"scale=5;\${STLJ[\$ave]} + (\${EARRAYLJ[\$run,\$ave]}*\${EARRAYLJ[\$run,\$ave]}) \" | bc -q)
AVQBND[\$ave]=\$( echo \"scale=5;\${AVQBND[\$ave]} + \${EARRAYQBND[\$run,\$ave]} \" | bc -q)
STQBND[\$ave]=\$( echo \"scale=5;\${STQBND[\$ave]} + (\${EARRAYQBND[\$run,\$ave]}*\${EARRAYQBND[\$run,\$ave]}) \" | bc -q)
AVQANG[\$ave]=\$( echo \"scale=5;\${AVQANG[\$ave]} + \${EARRAYQANG[\$run,\$ave]} \" | bc -q)
STQANG[\$ave]=\$( echo \"scale=5;\${STQANG[\$ave]} + (\${EARRAYQANG[\$run,\$ave]}*\${EARRAYQANG[\$run,\$ave]}) \" | bc -q)
AVQTOR[\$ave]=\$( echo \"scale=5;\${AVQTOR[\$ave]} + \${EARRAYQTOR[\$run,\$ave]} \" | bc -q)
STQTOR[\$ave]=\$( echo \"scale=5;\${STQTOR[\$ave]} + (\${EARRAYQTOR[\$run,\$ave]}*\${EARRAYQTOR[\$run,\$ave]}) \" | bc -q)
AVQIMP[\$ave]=\$( echo \"scale=5;\${AVQIMP[\$ave]} + \${EARRAYQIMP[\$run,\$ave]} \" | bc -q)
STQIMP[\$ave]=\$( echo \"scale=5;\${STQIMP[\$ave]} + (\${EARRAYQIMP[\$run,\$ave]}*\${EARRAYQIMP[\$run,\$ave]}) \" | bc -q)
AVBND[\$ave]=\$( echo \"scale=5;\${AVBND[\$ave]} + \${EARRAYBND[\$run,\$ave]} \" | bc -q)
STBND[\$ave]=\$( echo \"scale=5;\${STBND[\$ave]} + (\${EARRAYBND[\$run,\$ave]}*\${EARRAYBND[\$run,\$ave]}) \" | bc -q)
AVANG[\$ave]=\$( echo \"scale=5;\${AVANG[\$ave]} + \${EARRAYANG[\$run,\$ave]} \" | bc -q)
STANG[\$ave]=\$( echo \"scale=5;\${STANG[\$ave]} + (\${EARRAYANG[\$run,\$ave]}*\${EARRAYANG[\$run,\$ave]}) \" | bc -q)
AVTOR[\$ave]=\$( echo \"scale=5;\${AVTOR[\$ave]} + \${EARRAYTOR[\$run,\$ave]} \" | bc -q)
STTOR[\$ave]=\$( echo \"scale=5;\${STTOR[\$ave]} + (\${EARRAYTOR[\$run,\$ave]}*\${EARRAYTOR[\$run,\$ave]}) \" | bc -q)
AVIMP[\$ave]=\$( echo \"scale=5;\${AVIMP[\$ave]} + \${EARRAYIMP[\$run,\$ave]} \" | bc -q)
STIMP[\$ave]=\$( echo \"scale=5;\${STIMP[\$ave]} + (\${EARRAYIMP[\$run,\$ave]}*\${EARRAYIMP[\$run,\$ave]}) \" | bc -q)
AVTOT[\$ave]=\$( echo \"scale=5;\${AVTOT[\$ave]} + \${EARRAYTOT[\$run,\$ave]} \" | bc -q)
STTOT[\$ave]=\$( echo \"scale=5;\${STTOT[\$ave]} + (\${EARRAYTOT[\$run,\$ave]}*\${EARRAYTOT[\$run,\$ave]}) \" | bc -q)
AVKIN[\$ave]=\$( echo \"scale=5;\${AVKIN[\$ave]} + \${EARRAYKIN[\$run,\$ave]} \" | bc -q)
STKIN[\$ave]=\$( echo \"scale=5;\${STKIN[\$ave]} + (\${EARRAYKIN[\$run,\$ave]}*\${EARRAYKIN[\$run,\$ave]}) \" | bc -q)
AVPOT[\$ave]=\$( echo \"scale=5;\${AVPOT[\$ave]} + \${EARRAYPOT[\$run,\$ave]} \" | bc -q)
STPOT[\$ave]=\$( echo \"scale=5;\${STPOT[\$ave]} + (\${EARRAYPOT[\$run,\$ave]}*\${EARRAYPOT[\$run,\$ave]}) \" | bc -q)

done
done
echo \"#	QEL QVdW EL VdW QBND QANG QTOR QIMP BND ANG TOR IMP TOT KIN POT
# step  Lower bound     Upper bound  always\">$thisdir/${i}_${j}_${g}_benchmark.en
for ((ave=0;ave<\$LENGTH;ave++)) ; do

AVQEL[\$ave]=\$( echo \"scale=5;\${AVQEL[\$ave]} / $BENCHMARKS\" | bc -q)
AVQLJ[\$ave]=\$( echo \"scale=5;\${AVQLJ[\$ave]} / $BENCHMARKS\" | bc -q)
AVEL[\$ave]=\$( echo \"scale=5;\${AVEL[\$ave]} / $BENCHMARKS\" | bc -q)
AVLJ[\$ave]=\$( echo \"scale=5;\${AVLJ[\$ave]} / $BENCHMARKS\" | bc -q)
AVQBND[\$ave]=\$( echo \"scale=5;\${AVQBND[\$ave]} / $BENCHMARKS\" | bc -q)
AVQANG[\$ave]=\$( echo \"scale=5;\${AVQANG[\$ave]} / $BENCHMARKS\" | bc -q)
AVQTOR[\$ave]=\$( echo \"scale=5;\${AVQTOR[\$ave]} / $BENCHMARKS\" | bc -q)
AVQIMP[\$ave]=\$( echo \"scale=5;\${AVQIMP[\$ave]} / $BENCHMARKS\" | bc -q)
AVBND[\$ave]=\$( echo \"scale=5;\${AVBND[\$ave]} / $BENCHMARKS\" | bc -q)
AVANG[\$ave]=\$( echo \"scale=5;\${AVANG[\$ave]} / $BENCHMARKS\" | bc -q)
AVTOR[\$ave]=\$( echo \"scale=5;\${AVTOR[\$ave]} / $BENCHMARKS\" | bc -q)
AVIMP[\$ave]=\$( echo \"scale=5;\${AVIMP[\$ave]} / $BENCHMARKS\" | bc -q)
AVTOT[\$ave]=\$( echo \"scale=5;\${AVTOT[\$ave]} / $BENCHMARKS\" | bc -q)
AVKIN[\$ave]=\$( echo \"scale=5;\${AVKIN[\$ave]} / $BENCHMARKS\" | bc -q)
AVPOT[\$ave]=\$( echo \"scale=5;\${AVPOT[\$ave]} / $BENCHMARKS\" | bc -q)

STQEL[\$ave]=\$( echo \"scale=5;sqrt( (\${STQEL[\$ave]} / $BENCHMARKS) - ( \${AVQEL[\$ave]} * \${AVQEL[\$ave]} ))\" | bc -q)
STQLJ[\$ave]=\$( echo \"scale=5;sqrt( (\${STQLJ[\$ave]} / $BENCHMARKS) - ( \${AVQLJ[\$ave]} * \${AVQLJ[\$ave]} ))\" | bc -q)
STEL[\$ave]=\$( echo \"scale=5;sqrt( (\${STEL[\$ave]} / $BENCHMARKS) - ( \${AVEL[\$ave]} * \${AVEL[\$ave]} ))\" | bc -q)
STLJ[\$ave]=\$( echo \"scale=5;sqrt( (\${STLJ[\$ave]} / $BENCHMARKS) - ( \${AVLJ[\$ave]} * \${AVLJ[\$ave]} ))\" | bc -q)
STQBND[\$ave]=\$( echo \"scale=5;sqrt( (\${STQBND[\$ave]} / $BENCHMARKS) - ( \${AVQBND[\$ave]} * \${AVQBND[\$ave]} ))\" | bc -q)
STGANG[\$ave]=\$( echo \"scale=5;sqrt( (\${STQANG[\$ave]} / $BENCHMARKS) - ( \${AVQANG[\$ave]} * \${AVQANG[\$ave]} ))\" | bc -q)
STQTOR[\$ave]=\$( echo \"scale=5;sqrt( (\${STQTOR[\$ave]} / $BENCHMARKS) - ( \${AVQTOR[\$ave]} * \${AVQTOR[\$ave]} ))\" | bc -q)
STQIMP[\$ave]=\$( echo \"scale=5;sqrt( (\${STQIMP[\$ave]} / $BENCHMARKS) - ( \${AVQIMP[\$ave]} * \${AVQIMP[\$ave]} ))\" | bc -q)
STBND[\$ave]=\$( echo \"scale=5;sqrt( (\${STBND[\$ave]} / $BENCHMARKS) - ( \${AVBND[\$ave]} * \${AVBND[\$ave]} ))\" | bc -q)
STANG[\$ave]=\$( echo \"scale=5;sqrt( (\${STANG[\$ave]} / $BENCHMARKS) - ( \${AVANG[\$ave]} * \${AVANG[\$ave]} ))\" | bc -q)
STTOR[\$ave]=\$( echo \"scale=5;sqrt( (\${STTOR[\$ave]} / $BENCHMARKS) - ( \${AVTOR[\$ave]} * \${AVTOR[\$ave]} ))\" | bc -q)
STIMP[\$ave]=\$( echo \"scale=5;sqrt( (\${STIMP[\$ave]} / $BENCHMARKS) - ( \${AVIMP[\$ave]} * \${AVIMP[\$ave]} ))\" | bc -q)
STTOT[\$ave]=\$( echo \"scale=5;sqrt( (\${STTOT[\$ave]} / $BENCHMARKS) - ( \${AVTOT[\$ave]} * \${AVTOT[\$ave]} ))\" | bc -q)
STKIN[\$ave]=\$( echo \"scale=5;sqrt( (\${STKIN[\$ave]} / $BENCHMARKS) - ( \${AVKIN[\$ave]} * \${AVKIN[\$ave]} ))\" | bc -q)
STPOT[\$ave]=\$( echo \"scale=5;sqrt( (\${STPOT[\$ave]} / $BENCHMARKS) - ( \${AVPOT[\$ave]} * \${AVPOT[\$ave]} ))\" | bc -q)


UPPERQEL[\$ave]=\$( echo \"scale=3;\${AVQEL[\$ave]} + \${STQEL[\$ave]}\" | bc -q)
LOWERQEL[\$ave]=\$( echo \"scale=3;\${AVQEL[\$ave]} - \${STQEL[\$ave]}\" | bc -q)
UPPERQLJ[\$ave]=\$( echo \"scale=3;\${AVQLJ[\$ave]} + \${STQLJ[\$ave]}\" | bc -q)
LOWERQLJ[\$ave]=\$( echo \"scale=3;\${AVQLJ[\$ave]} - \${STQLJ[\$ave]}\" | bc -q)
UPPEREL[\$ave]=\$( echo \"scale=3;\${AVEL[\$ave]} + \${STEL[\$ave]}\" | bc -q)
LOWEREL[\$ave]=\$( echo \"scale=3;\${AVEL[\$ave]} - \${STEL[\$ave]}\" | bc -q)
UPPERLJ[\$ave]=\$( echo \"scale=3;\${AVLJ[\$ave]} + \${STLJ[\$ave]}\" | bc -q)
LOWERLJ[\$ave]=\$( echo \"scale=3;\${AVLJ[\$ave]} - \${STLJ[\$ave]}\" | bc -q)
UPPERQBND[\$ave]=\$( echo \"scale=3;\${AVQBND[\$ave]} + \${STQBND[\$ave]}\" | bc -q)
LOWERQBND[\$ave]=\$( echo \"scale=3;\${AVQBND[\$ave]} - \${STQBND[\$ave]}\" | bc -q)
UPPERQANG[\$ave]=\$( echo \"scale=3;\${AVQANG[\$ave]} + \${STQANG[\$ave]}\" | bc -q)
LOWERQANG[\$ave]=\$( echo \"scale=3;\${AVQANG[\$ave]} - \${STQANG[\$ave]}\" | bc -q)
UPPERQTOR[\$ave]=\$( echo \"scale=3;\${AVQTOR[\$ave]} + \${STQTOR[\$ave]}\" | bc -q)
LOWERQTOR[\$ave]=\$( echo \"scale=3;\${AVQTOR[\$ave]} - \${STQTOR[\$ave]}\" | bc -q)
UPPERQIMP[\$ave]=\$( echo \"scale=3;\${AVQIMP[\$ave]} + \${STQIMP[\$ave]}\" | bc -q)
LOWERQIMP[\$ave]=\$( echo \"scale=3;\${AVQIMP[\$ave]} - \${STQIMP[\$ave]}\" | bc -q)
UPPERBND[\$ave]=\$( echo \"scale=3;\${AVBND[\$ave]} + \${STBND[\$ave]}\" | bc -q)
LOWERBND[\$ave]=\$( echo \"scale=3;\${AVBND[\$ave]} - \${STBND[\$ave]}\" | bc -q)
UPPERANG[\$ave]=\$( echo \"scale=3;\${AVANG[\$ave]} + \${STANG[\$ave]}\" | bc -q)
LOWERANG[\$ave]=\$( echo \"scale=3;\${AVANG[\$ave]} - \${STANG[\$ave]}\" | bc -q)
UPPERTOR[\$ave]=\$( echo \"scale=3;\${AVTOR[\$ave]} + \${STTOT[\$ave]}\" | bc -q)
LOWERTOR[\$ave]=\$( echo \"scale=3;\${AVTOR[\$ave]} - \${STTOT[\$ave]}\" | bc -q)
UPPERIMP[\$ave]=\$( echo \"scale=3;\${AVIMP[\$ave]} + \${STIMP[\$ave]}\" | bc -q)
LOWERIMP[\$ave]=\$( echo \"scale=3;\${AVIMP[\$ave]} - \${STIMP[\$ave]}\" | bc -q)
UPPERTOT[\$ave]=\$( echo \"scale=3;\${AVTOT[\$ave]} + \${STTOT[\$ave]}\" | bc -q)
LOWERTOT[\$ave]=\$( echo \"scale=3;\${AVTOT[\$ave]} - \${STTOT[\$ave]}\" | bc -q)
UPPERKIN[\$ave]=\$( echo \"scale=3;\${AVKIN[\$ave]} + \${STKIN[\$ave]}\" | bc -q)
LOWERKIN[\$ave]=\$( echo \"scale=3;\${AVKIN[\$ave]} - \${STKIN[\$ave]}\" | bc -q)
UPPERPOT[\$ave]=\$( echo \"scale=3;\${AVPOT[\$ave]} + \${STPOT[\$ave]}\" | bc -q)
LOWERPOT[\$ave]=\$( echo \"scale=3;\${AVPOT[\$ave]} - \${STPOT[\$ave]}\" | bc -q)

NUM=\$(( \$ave + 1 ))
echo \"\$NUM	\${LOWERQEL[\$ave]} \${UPPERQEL[\$ave]} \${LOWERQLJ[\$ave]} \${UPPERQLJ[\$ave]} \${LOWEREL[\$ave]} \${UPPEREL[\$ave]} \${LOWERLJ[\$ave]} \${UPPERLJ[\$ave]} \${LOWERQBND[\$ave]} \${UPPERQBND[\$ave]} \${LOWERQANG[\$ave]} \${UPPERQANG[\$ave]} \${LOWERQTOR[\$ave]} \${UPPERQTOR[\$ave]} \${LOWERQIMP[\$ave]} \${UPPERQIMP[\$ave]} \${LOWERBND[\$ave]} \${UPPERBND[\$ave]} \${LOWERANG[\$ave]} \${UPPERANG[\$ave]} \${LOWERTOR[\$ave]} \${UPPERTOR[\$ave]} \${LOWERIMP[\$ave]} \${UPPERIMP[\$ave]} \${LOWERTOT[\$ave]} \${UPPERTOT[\$ave]} \${LOWERKIN[\$ave]} \${UPPERKIN[\$ave]} \${LOWERPOT[\$ave]} \${UPPERPOT[\$ave]}\">>$thisdir/${i}_${j}_${g}_benchmark.en
done
">>basetest.sh
$subcommand basetest.sh
cd $thisdir
fi
done
done
done


}

#now comes the part where the test function is called with the right binary as argument
make_dirs $1


