#!/bin/bash 
#################################################################
# NOTE:
# Uncomment or modify the next two lines depending on your system
# openmpi or impi version, and the path to your Q binaries.
#################################################################
# module load openmpi-x86_64
# export bindir="/home/username/software/qsource/bin"

if [ "x$bindir" == "x" ]
then 
 echo "Please set the bindir variable to point to the Q directory."
 exit 1
elif [ ! -x $bindir/qdynp ]
then
 echo "Can't locate qdynp in the bindir, or you don't have executable permission"
 exit 1
else
 echo "Detected qdynp in ${bindir}"
fi

# How many cores on this machine?
# grep "cpu cores" /proc/cpuinfo  
# cpu cores: 4
# The cores need then to be added to get the total number of cores available per node.
# For now bc is doing the sum, but BEWARE, maybe bc is not installed on all nodes.
# CORES=`grep processor /proc/cpuinfo | wc -l`
# CORES=`grep cores /proc/cpuinfo | awk '{print $4}' | paste -sd+ | bc`
CORES=8
echo "Running simulation on $CORES cores."

rm -f eq{1..5}.log dc{1..5}.log >& /dev/null

# Useful vars
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

for step in {1..5}
do
 echo -n "Running equilibration step ${step} of 5                         "
 if mpirun -np $CORES $bindir/qdynp eq${step}.inp > eq${step}.log
 then echo -e "$OK"
 else 
  echo -e "$FAILED"
  echo "Check output (eq${step}.log) for more info."
  exit 1
 fi
done


for step in {1..5}
do
 echo -n "Running production run step ${step} of 5                        "
 if mpirun -np $CORES $bindir/qdynp dc${step}.inp > dc${step}.log
  then echo -e "$OK"
 else 
  echo -e "$FAILED"
  echo "Check output (dc${step}.logl) for more info."
  exit 1
 fi
done


