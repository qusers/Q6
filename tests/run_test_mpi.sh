#!/bin/bash 

if [ "x$Q_PATH" == "x" ]
then 
 echo "Please set the Q_PATH variable to point at the Q directory."
 exit 1
elif [ ! -x $Q_PATH/Qdyn5p ]
then
 echo "Can't locate Qdyn5P in the Q_PATH, or you don't have 
       execute permisson."
 exit 1
else
 echo "Detected Qdyn5p in ${Q_PATH}"
fi

# How many cores on this machine?
# Ex:
# grep processor /proc/cpuinfo  
#	processor	: 0
#	processor	: 1
#	processor	: 2
#	processor	: 3
# Numer of lines is 4, so four cores in this case.
CORES=`grep processor /proc/cpuinfo | wc -l`
echo "Running simulation on $CORES cores."

rm eq{1..5}.log dc{1..5}.log >& /dev/null

# Useful vars
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

for step in {1..5}
do
 echo -n "Running equilibration step ${step} of 5                         "
 if mpirun -n $CORES $Q_PATH/Qdyn5p eq${step}.inp > eq${step}.log
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
 if mpirun -n $CORES $Q_PATH/Qdyn5p dc${step}.inp > dc${step}.log
  then echo -e "$OK"
 else 
  echo -e "$FAILED"
  echo "Check output (dc${step}.log) for more info."
  exit 1
 fi
done


