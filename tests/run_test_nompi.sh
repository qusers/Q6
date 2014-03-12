#!/bin/bash 

qbinary=$QDIR/bin/qdyn5

if [ "x$QDIR" == "x" ]
then 
 echo "Please set the QDIR variable to point to the Q folder"
 exit 1
elif [ ! -x $qbinary ]
then
 echo "Can't locate qdyn in the QDIR, or you don't have 
       execute permisson."
 exit 1
else
 echo "Detected qdyn in ${QDIR}"
fi

rm eq{1..5}.log dc{1..5}.log >& /dev/null

# Useful vars
OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

for step in {1..5}
do
 echo -n "Running equilibration step ${step} of 5                         "
 if $qbinary eq${step}.inp > eq${step}.log
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
 if $qbinary dc${step}.inp > dc${step}.log
  then echo -e "$OK"
 else 
  echo -e "$FAILED"
  echo "Check output (dc${step}.log) for more info."
  exit 1
 fi
done


