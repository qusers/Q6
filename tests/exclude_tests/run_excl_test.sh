#!/bin/bash

#script is for running locally, for cluster please add #SBATCH information
#takes arguments parallel serial hybrid relax
#runs both system with exclusion and without for parallel or serial
#relax creates the restart file again, will take a long time
#uses parallel version for that

#written by paul.bauer@icm.uu.se
#released under Beer-Ware Licence

QPATH="/data/work/q_source/bin"
CORES=4
OMPCORES=2

STARTDIR=`pwd`

MPICOM="mpirun "
if [ $# -ne 1 ] ; then
echo "Number of arguments has to be one for this script to work"
echo "Please supply what test should be run"
echo "Supported arguments are hybrid  parallel serial relax"
exit 1
fi
if [ "$1" != "parallel" ] && [ "$1" != "serial" ] && [ "$1" != "relax" ]  && [ "$1" != "hybrid" ]
then
echo "Did not understand the argument $1"
echo "Please supply one of the following: serial parallel hybrid relax"
exit 2
fi

#This makes the actual run test script part
#arguments are file name and type of run
function make_run_command() {
if [ "$1" == "serial" ]  ; then
QBIN=$QPATH/qdyn5
RUNCOMA="$QBIN"
RUNCOMB=""
fi
if [ "$1" == "parallel" ] || [ "$1" == "relax" ] ; then
QBIN=$QPATH/qdyn5p
RUNCOMA="$MPICOM -np $CORES $QBIN"
RUNCOMB=""
fi
if [ "$1" == "hybrid" ] ; then
QBIN=$QPATH/qdyn5h
RUNCOMA="$MPICOM -np 1 -x OMP_NUM_THREADS=2 $QBIN " 
RUNCOMB=": -np \$(( $CORES - 2 )) -x OMP_NUM_THREADS=1 $QBIN"
fi


}

#function that runs test EVB MAPPING trajectory
function run_test() {

make_run_command $1

if [ "$1" != "relax" ] 
then
for kind in excl noexcl 
do

#removing oldes files
#keeping one backup
if [ -d run_${1}_${kind} ]
then
  if [ -d run_${1}_${kind}_bkp ]
  then
    echo "removing old ${1}_${kind} backup"
    /bin/rm -rf run_${1}_${kind}_bkp
  fi
  echo "backing up old ${1}_${kind}"
  mv run_${1}_${kind} run_${1}_${kind}_bkp
fi


mkdir run_${1}_${kind}
cp $STARTDIR/inputs/${kind}/*pl $STARTDIR/inputs/*fep $STARTDIR/inputs/*top $STARTDIR/inputs/*re run_${1}_${kind}
cd run_${1}_${kind}

if [[ $SNIC_TMP != "" ]];then
MYTMPDIR=$SNIC_TMP/run
else
#Most likely a local run
MYTMPDIR=$$
MYTMPDIR=/dev/shm/${MYTMPDIR}_local
fi
THISDIR=${PWD}

perl gen_inps.pl
for (( c=100; c>=0; c-- )); do
NUM=${c}
if [[ ${c} -lt 100 ]];then
NUM=0${c}
fi
if [[ ${c} -lt 10 ]];then
NUM=00${c}
fi
MODULUS=$(( ${c}%2 ))
if [[ $MODULUS = 0 ]];then
mkdir ${MYTMPDIR}
 echo "Running FEP EVB Mapping frame for ${NUM}"
 cp *top *re *fep fep_${NUM}.inp ${MYTMPDIR}
 cd ${MYTMPDIR}
 echo "Entering temporary directory"
 $RUNCOMA fep_${NUM}.inp >fep_${NUM}.log $RUNCOMB
 if [ `grep "Terminated normally" fep_${NUM}.log | wc -l` -ne 0 ]
 then echo -e "$OK\nCopying back the files"
 cp * ${THISDIR}
 cd ${THISDIR}
 rm -r ${MYTMPDIR}
 else
  echo -e "$FAILED"
 echo "Check output (fep_${NUM}.log) for more info."
  echo "Copying back the files"
  cp * ${THISDIR}
  cd ${THISDIR}
  rm -r ${MYTMPDIR}
  exit 1
 fi
fi
done

cd $STARTDIR
done


else
#removing oldes files
#keeping one backup
if [ -d run_relax ]
then
  if [ -d run_relax_bkp ] 
  then
    echo "removing old relax backup"
    /bin/rm -rf run_relax_bkp
  fi
  echo "backing up old relax"
  mv run_relax run_relax_bkp
fi


mkdir run_relax
cp inputs/relax/*inp inputs/relax/*top inputs/relax/*fep run_relax
cd run_relax
if [[ $SNIC_TMP != "" ]];then
MYTMPDIR=$SNIC_TMP/run
else
#Most likely a local run
MYTMPDIR=$$
MYTMPDIR=/dev/shm/${MYTMPDIR}_local
fi
THISDIR=${PWD}
for step in  1-dyn_rlx_wat_noshake_verysmall 2-dyn_rlx_wat_noshake_bitsmall 3-dyn_rlx_wat_noshake_small 4-dyn_rlx_wat_noshake 5-dyn_rlx_wat 6-dyn_rlx_wat_eq 7-dyn_cool 8-dyn_RLS 9-dyn_warm150k 10-dyn_warmRT 11-dyn_eqRT 12-dyn_eqPUB
do
 mkdir ${MYTMPDIR}
 echo "Running equilibration step ${step}"
 cp ${step}.inp *fep *top *re ${MYTMPDIR}
 cd ${MYTMPDIR}
 echo "Entering temporary directory"
 echo "Now in ${MYTMPDIR}"
 $RUNCOMA  ${step}.inp > ${step}.log $RUNCOMB
 if [ `grep "Terminated normally" ${step}.log | wc -l` -ne 0 ]
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

cd $STARTDIR
fi
}


run_test $1


