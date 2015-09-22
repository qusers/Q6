#!/bin/bash -l
#SBATCH -J relax
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -A SNIC2013-26-1

#module add intel/12.1.4 impi/4.0.3.008

Q_PATH="/home/acmnpv/data/work/q_source/bin"
if [[ $SNIC_TMP != "" ]];then
MYTMPDIR=$SNIC_TMP/run


else
#Most likely a local run
MYTMPDIR=$$
MYTMPDIR=/dev/shm/${MYTMPDIR}_local

fi

#mkdir ${MYTMPDIR}
THISDIR=${PWD}

THISCORES=2

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"



if [ "x$Q_PATH" == "x" ]
then
 echo "Please set the Q_PATH variable to point at the Q directory."
 exit 1
elif [ ! -x $Q_PATH/qdyn5p ]
then
 echo "Can't locate qdyn5p in the Q_PATH, or you don't have 
       execute permisson."
 exit 1
else
 echo "Detected qdyn5p in ${Q_PATH}"
fi

echo "Running simulation on $THISCORES cores."

#cp ../epox_int_lig.fep lig.fep
#cp ../*top .
for step in  1-dyn_rlx_wat_noshake_verysmall 2-dyn_rlx_wat_noshake_bitsmall 3-dyn_rlx_wat_noshake_small 4-dyn_rlx_wat_noshake 5-dyn_rlx_wat 6-dyn_rlx_wat_eq 7-dyn_cool 8-dyn_RLS 9-dyn_warm150k 10-dyn_warmRT 11-dyn_eqRT 12-dyn_eqPUB
do
 mkdir ${MYTMPDIR}
 echo "Running equilibration step ${step}"
 cp ${step}.inp *fep *top *re ${MYTMPDIR}
 cd ${MYTMPDIR}
 echo "Entering temporary directory"
 echo "Now in ${MYTMPDIR}"
 if mpirun --bind-to-core --bycore --report-bindings -n ${THISCORES} $Q_PATH/qdyn5p ${step}.inp > ${step}.log
 then echo -e "$OK\nCopying back the files"
 cp * ${THISDIR}
 cd ${THISDIR}
 rm -r ${MYTMPDIR}
 cp ${step}.re ${step}_rest.re
 else
  echo -e "$FAILED"
  echo "Check output (${step}.log) for more info."
  echo "Copying back the files"
  cp * ${THISDIR}
  cd ${THISDIR}
  rm -r ${MYTMPDIR}
  exit 1
 fi
done

#cd ..
#cp relaxation_it/12-dyn_eqPUB.re relaxation_it/2cjpFH_ionres_oplsa.top relaxation_it/lig.fep framegen_it

#cd framegen_it
#cp 12-dyn_eqPUB.re 0.re
#cp 0.re 0_rest.re

#sbatch frame_gen.sh
