#!/bin/bash -l
#SBATCH -J eq_water
#SBATCH -p devel -n 8
#SBATCH -t 00:60:00 
#SBATCH -A s00212-26 

cp ../*top .

module load intel openmpi
mpirun Qdyn5 fep_warm20k.inp>fep_warm20k.log
mpirun Qdyn5 fep_warm100k.inp>fep_warm100k.log
mpirun Qdyn5 fep_warm300k.inp>fep_warm300k.log


mpirun Qdyn5 fep_100.inp >fep_100.log
mpirun Qdyn5 fep_098.inp >fep_098.log
mpirun Qdyn5 fep_096.inp >fep_096.log
mpirun Qdyn5 fep_094.inp >fep_094.log
mpirun Qdyn5 fep_092.inp >fep_092.log
mpirun Qdyn5 fep_090.inp >fep_090.log
mpirun Qdyn5 fep_088.inp >fep_088.log
mpirun Qdyn5 fep_086.inp >fep_086.log
mpirun Qdyn5 fep_084.inp >fep_084.log
mpirun Qdyn5 fep_082.inp >fep_082.log
mpirun Qdyn5 fep_080.inp >fep_080.log
mpirun Qdyn5 fep_078.inp >fep_078.log
mpirun Qdyn5 fep_076.inp >fep_076.log
mpirun Qdyn5 fep_074.inp >fep_074.log
mpirun Qdyn5 fep_072.inp >fep_072.log
mpirun Qdyn5 fep_070.inp >fep_070.log
mpirun Qdyn5 fep_068.inp >fep_068.log
mpirun Qdyn5 fep_066.inp >fep_066.log
mpirun Qdyn5 fep_064.inp >fep_064.log
mpirun Qdyn5 fep_062.inp >fep_062.log
mpirun Qdyn5 fep_060.inp >fep_060.log
mpirun Qdyn5 fep_058.inp >fep_058.log
mpirun Qdyn5 fep_056.inp >fep_056.log
mpirun Qdyn5 fep_054.inp >fep_054.log
mpirun Qdyn5 fep_052.inp >fep_052.log
mpirun Qdyn5 fep_050.inp >fep_050.log
mpirun Qdyn5 fep_048.inp >fep_048.log
mpirun Qdyn5 fep_046.inp >fep_046.log
mpirun Qdyn5 fep_044.inp >fep_044.log
mpirun Qdyn5 fep_042.inp >fep_042.log
mpirun Qdyn5 fep_040.inp >fep_040.log
mpirun Qdyn5 fep_038.inp >fep_038.log
mpirun Qdyn5 fep_036.inp >fep_036.log
mpirun Qdyn5 fep_034.inp >fep_034.log
mpirun Qdyn5 fep_032.inp >fep_032.log
mpirun Qdyn5 fep_030.inp >fep_030.log
mpirun Qdyn5 fep_028.inp >fep_028.log
mpirun Qdyn5 fep_026.inp >fep_026.log
mpirun Qdyn5 fep_024.inp >fep_024.log
mpirun Qdyn5 fep_022.inp >fep_022.log
mpirun Qdyn5 fep_020.inp >fep_020.log
mpirun Qdyn5 fep_018.inp >fep_018.log
mpirun Qdyn5 fep_016.inp >fep_016.log
mpirun Qdyn5 fep_014.inp >fep_014.log
mpirun Qdyn5 fep_012.inp >fep_012.log
mpirun Qdyn5 fep_010.inp >fep_010.log
mpirun Qdyn5 fep_008.inp >fep_008.log
mpirun Qdyn5 fep_006.inp >fep_006.log
mpirun Qdyn5 fep_004.inp >fep_004.log
mpirun Qdyn5 fep_002.inp >fep_002.log
mpirun Qdyn5 fep_000.inp >fep_000.log
