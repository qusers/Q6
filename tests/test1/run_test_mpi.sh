#!/bin/bash -l
#SBATCH -J conser
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH -A SNIC2013-26-1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=paul.bauer@icm.uu.se
#################################################################
# NOTE:
# Uncomment or modify the next two lines depending on your system
# openmpi or impi version, and the path to your Q binaries.
#################################################################


module load intel/12.1.4 impi/4.0.3.008

#set -e
QDIR=/home/x_pauba/glob/debug/noseho/bin

qbinary=$QDIR/qdyn5p
wd=`pwd`

#set -e

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

function run_test() {
rm -f eq{1..5}.log dc{1..6}.log >& /dev/null
tdir="${MYTMPDIR}_$1"
cdir=`pwd`
for step in {1..5}
do
 mkdir ${tdir}
 echo "Running equilibration step ${step}"
 cp eq${step}.inp *fep *top *re ${tdir}
 cd ${tdir}
 echo "Entering temporary directory"
 echo "Now in ${tdir}"
 mpprun --nranks=$CORES $qbinary eq${step}.inp > eq${step}.log
 check=$( grep "terminated normally" eq${step}.log | wc -l)
 if [ ! $check -eq 0 ]
 then echo -e "$OK\nCopying back the files"
 cp * ${cdir}
 cd ${cdir}
 rm -r ${tdir}
 else
  echo -e "$FAILED"
  echo "Check output (eq${step}.log) for more info."
  echo "Copying back the files"
  cp * ${cdir}
  cd ${cdir}
  rm -r ${tdir}
  exit 1
 fi
done

for step in {1..6}
do
 mkdir ${tdir}
 echo "Running production run step ${step} of 6"
 cp dc${step}.inp *fep *top *re ${tdir}
 cd ${tdir}
 echo "Entering temporary directory"
 echo "Now in ${tdir}"
 mpprun --nranks=$CORES $qbinary dc${step}.inp > dc${step}.log
 check=$( grep "terminated normally" dc${step}.log | wc -l)
 if [ ! $check -eq 0 ]
 then echo -e "$OK\nCopying back the files"
 cp * ${cdir}
 cd ${cdir}
 rm -r ${tdir}
 else
  echo -e "$FAILED"
  echo "Check output (dc${step}.log) for more info"
  echo "Copying back the files"
  cp * ${cdir}
  cd ${cdir}
  rm -r ${tdir}
  exit 1
 fi
done
}


# How many cores on this machine?
#  grep "cpu cores" /proc/cpuinfo  
# cpu cores: 4
# cpu cores: 4
# cpu cores: 4
# cpu cores: 4
# The cores need then to be added to get the total number of cores available per node.
# For now bc is doing the sum, but BEWARE, maybe bc is not installed in all
# nodes.
#CORES=`grep processor /proc/cpuinfo | wc -l`
CORES=$SLURM_NPROCS
echo "Running simulation on $CORES cores."

rm -rf $wd/run-mpitest
mkdir $wd/run-mpitest

cp *inp eval_test.sh qsurr_benchmark.en $wd/run-mpitest
cp $wd/prep/lig_w.top $wd/run-mpitest/lig_w.top
cp $wd/prep/lig_w.fep $wd/run-mpitest/lig_w.fep

cd $wd/run-mpitest

run_test mpi

cd ..
