#!/bin/bash  -l
#SBATCH -J SPH_LAN_benvv
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -A SNIC2013-26-1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=paul.bauer@icm.uu.se
#make benchmark
#run 20 different random seeds to get some variation
#by default, 4 tests are run at the same time
#using the inbuild queue system
#paul.bauer@icm.uu.se
#as always, beer-ware license

module load intel/12.1.4 impi/4.0.3.008

#set -e
QDIR=/home/x_pauba/glob/debug/noseho/bin
BENCHMARKS=20
PARALLEL=14
RUNNING=0
wd=`pwd`
qbinary=$QDIR/qdyn5
# Useful vars
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

if [[ $SNIC_TMP != "" ]];then
MYTMPDIR=$SNIC_TMP/run_$$


else
#Most likely a local run
MYTMPDIR=$$
MYTMPDIR=/dev/shm/${MYTMPDIR}_local

fi


if [ -z $qbinary ]
then
 echo "Please set the qbinary variable to point to the Q folder"
 exit 1
elif [ ! -x $qbinary ]
then
 echo "Can't locate qdyn in the $qbinary variable, or you don't have
       execute permisson."
 exit 1
else
 echo "Detected qdyn in ${QDIR}"
fi


function run_test() {
rm -f eq{1..5}.log dc{1..6}.log >& /dev/null
tdir[$1]="${MYTMPDIR}_$1-benchmark"
cdir[$1]=`pwd`
for step in {1..5}
do
 mkdir ${tdir[$1]}
 echo "Running equilibration step ${step}"
 cp eq${step}.inp *fep *top *re ${tdir[$1]}
 cd ${tdir[$1]}
 echo "Entering temporary directory"
 echo "Now in ${tdir[$1]}"
 $qbinary eq${step}.inp > eq${step}.log
 check=$( grep "terminated normally" eq${step}.log | wc -l)
 if [ ! $check -eq 0 ]
 then echo -e "$OK\nCopying back the files"
 cp * ${cdir[$1]}
 cd ${cdir[$1]}
 rm -r ${tdir[$1]}
 else
  echo -e "$FAILED"
  echo "Check output (eq${step}.log) for more info."
  echo "Copying back the files"
  cp * ${cdir[$1]}
  cd ${cdir[$1]}
  rm -r ${tdir[$1]}
  exit 1
 fi
done

for step in {1..6}
do
 mkdir ${tdir[$1]}
 echo "Running production run step ${step} of 6"
 cp dc${step}.inp *fep *top *re ${tdir[$1]}
 cd ${tdir[$1]}
 echo "Entering temporary directory"
 echo "Now in ${tdir[$1]}"
 $qbinary dc${step}.inp > dc${step}.log
 check=$( grep "terminated normally" dc${step}.log | wc -l)
 if [ ! $check -eq 0 ]
 then echo -e "$OK\nCopying back the files"
 cp * ${cdir[$1]}
 cd ${cdir[$1]}
 rm -r ${tdir[$1]}
 else
  echo -e "$FAILED"
  echo "Check output (dc${step}.log) for more info"
  echo "Copying back the files"
  cp * ${cdir[$1]}
  cd ${cdir[$1]}
  rm -r ${tdir[$1]}
  exit 1
 fi
done
}

for ((run=1;run<=$BENCHMARKS;run++)) ; do
rm -rf $wd/run-$run
mkdir $wd/run-$run
cp *inp $wd/run-$run
 
cp $wd/prep/lig_w.top $wd/run-$run/lig_w.top
cp $wd/prep/lig_w.fep $wd/run-$run/lig_w.fep

sed -i s/"random_seed.*1"/"random_seed                $RANDOM"/g $wd/run-$run/eq1.inp



cd $wd/run-$run

run_test $run &>logfile & 
RARRAY[$RUNNING]=$$
RUNNING=$(( ${RUNNING} + 1 ))
cd $wd

if [ $RUNNING -eq $PARALLEL ] || [ $run -eq $BENCHMARKS ] ; then
wait 
for ((testr=$run;testr>=$(($run - $RUNNING + 1));testr--)) ; do
	err=`grep "terminated normally" $wd/run-$testr/dc6.log | wc -l`
	if [ $err -eq 0 ] ; then
	echo "Error during running of $wd/run-$testr"
	echo "Aborting tests"
	exit 666
	fi
done

RUNNING=0
fi

done


#Now get energies in nice arrays to prepare the new benchmark file
#TARRAY=( `grep "Q-surr. 1 1.0000" $wd/run-1/eq*log $wd/run-1/dc*log | sed s/"Q-surr. 1 1.0000"//g | awk '{print $1}'` )
#LENGTH=${#TARRAY[@]}
LENGTH=`grep "Q-surr. 1 1.0000" $wd/run-1/eq*log $wd/run-1/dc*log | sed s/".*Q-surr. 1 1.0000"//g | awk '{print $1}' | wc -l`
unset AVEL STEL AVLJ STLJ
for ((ave=0;ave<$LENGTH;ave++)) ; do
AVEL[$ave]=0
STEL[$ave]=0
AVLJ[$ave]=0
STLJ[$ave]=0
done
for ((run=1;run<=$BENCHMARKS;run++)) ; do
unset TARRAYEL TARRAYLJ
TARRAYEL=( `grep "Q-surr.*1 1.0000" $wd/run-$run/eq*log $wd/run-$run/dc*log | sed s/".*Q-surr. 1 1.0000"//g | awk '{print $1}'` )
TARRAYLJ=( `grep "Q-surr.*1 1.0000" $wd/run-$run/eq*log $wd/run-$run/dc*log | sed s/".*Q-surr. 1 1.0000"//g | awk '{print $2}'` )
for ((ave=0;ave<$LENGTH;ave++)) ; do
EARRAYEL[$run,$ave]=${TARRAYEL[$ave]}
EARRAYLJ[$run,$ave]=${TARRAYLJ[$ave]}
done
for ((ave=0;ave<$LENGTH;ave++)) ; do
AVEL[$ave]=$( echo "scale=5;${AVEL[$ave]} + ${EARRAYEL[$run,$ave]} " | bc -q)
STEL[$ave]=$( echo "scale=5;${STEL[$ave]} + (${EARRAYEL[$run,$ave]}*${EARRAYEL[$run,$ave]}) " | bc -q)
AVLJ[$ave]=$( echo "scale=5;${AVLJ[$ave]} + ${EARRAYLJ[$run,$ave]} " | bc -q)
STLJ[$ave]=$( echo "scale=5;${STLJ[$ave]} + (${EARRAYLJ[$run,$ave]}*${EARRAYLJ[$run,$ave]}) " | bc -q)
done
done
echo "#       el              vdw     
# step  Lower bound     Upper bound     Lower bound     Upper bound">$wd/qsurr_benchmark.en
for ((ave=0;ave<$LENGTH;ave++)) ; do

AVEL[$ave]=$( echo "scale=5;${AVEL[$ave]} / $BENCHMARKS" | bc -q)
AVLJ[$ave]=$( echo "scale=5;${AVLJ[$ave]} / $BENCHMARKS" | bc -q)


STEL[$ave]=$( echo "scale=5;sqrt( (${STEL[$ave]} / $BENCHMARKS) - ( ${AVEL[$ave]} * ${AVEL[$ave]} ))" | bc -q)
STLJ[$ave]=$( echo "scale=5;sqrt( (${STLJ[$ave]} / $BENCHMARKS) - ( ${AVLJ[$ave]} * ${AVLJ[$ave]} ))" | bc -q) 
UPPEREL[$ave]=$( echo "scale=3;${AVEL[$ave]} + ${STEL[$ave]}" | bc -q)
LOWEREL[$ave]=$( echo "scale=3;${AVEL[$ave]} - ${STEL[$ave]}" | bc -q)
UPPERLJ[$ave]=$( echo "scale=3;${AVLJ[$ave]} + ${STLJ[$ave]}" | bc -q)
LOWERLJ[$ave]=$( echo "scale=3;${AVLJ[$ave]} - ${STLJ[$ave]}" | bc -q)
NUM=$(( $ave + 1 ))
echo "$NUM	${LOWEREL[$ave]} ${UPPEREL[$ave]} ${LOWERLJ[$ave]} ${UPPERLJ[$ave]}">>$wd/qsurr_benchmark.en
done


